/*
 * LI_Calib: An Open Platform for LiDAR-IMU Calibration
 * Copyright (C) 2020 Jiajun Lv
 * Copyright (C) 2020 Kewei Hu
 * Copyright (C) 2020 Jinhong Xu
 * Copyright (C) 2020 LI_Calib Contributors
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */
#include <core/surfel_association.h>
#include <pcl/common/common.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/extract_clusters.h>
#include <ctime>
#include <cmath>
#include "omp.h"


namespace licalib {

void SurfelAssociation::initColorList() {
  color_list_.set_capacity(6);
  color_list_.push_back(0xFF0000);
  color_list_.push_back(0xFF00FF);
  color_list_.push_back(0x436EEE);
  color_list_.push_back(0xBF3EFF);
  color_list_.push_back(0xB4EEB4);
  color_list_.push_back(0xFFE7BA);
}

void SurfelAssociation::clearSurfelMap() {
  surfel_planes_.clear();
  spoint_per_surfel_.clear();
  spoint_downsampled_.clear();
  surfels_map_.clear();
  spoints_all_.clear();
}

// 每个格子拟合平面，并且给不同cell的内点不同颜色
void SurfelAssociation::setSurfelMap(
        const pclomp::NormalDistributionsTransform<VPoint, VPoint>::Ptr& ndtPtr,
        double timestamp) {
  clearSurfelMap();
  map_timestamp_ = timestamp;

  //mCellSize = ndtPtr->getTargetCells().getLeafSize()(0);

  // check each cell
  Eigen::Vector3i counter(0,0,0);
  for (const auto &v : ndtPtr->getTargetCells().getLeaves()) { //leaves报保存的是map<int,leaf> leaf保存当前cell的信息
    auto leaf = v.second;

    if (leaf.nr_points < 10) // cell点数要>10
      continue;
    // leaf直接可以获取
    int plane_type = checkPlaneType(leaf.getEvals(), leaf.getEvecs(), p_lambda_); //cell里的特征值和特征向量，用于判断是否为平面
    if (plane_type < 0)
      continue;

    Eigen::Vector4d surfCoeff;
    VPointCloud::Ptr cloud_inliers = VPointCloud::Ptr(new VPointCloud);
    if (!fitPlane(leaf.pointList_.makeShared(), surfCoeff, cloud_inliers)) // cell里的点RANSAC拟合平面，平面参数保存到surfCoeff中，平面内点保存到cloud_inliers中，inlier点数必须大与20
      continue;

    counter(plane_type) += 1; // -1 0 1 2，与坐标轴平行的序号
    // 平面信息保存
    SurfelPlane surfplane;
    surfplane.cloud = leaf.pointList_;
    surfplane.cloud_inlier = *cloud_inliers;
    surfplane.p4 = surfCoeff;
    surfplane.Pi = -surfCoeff(3) * surfCoeff.head<3>();

    // 平面包围框
    VPoint min, max;
    pcl::getMinMax3D(surfplane.cloud, min, max);
    surfplane.boxMin = Eigen::Vector3d(min.x, min.y, min.z);
    surfplane.boxMax = Eigen::Vector3d(max.x, max.y, max.z);

    surfel_planes_.push_back(surfplane);
  }

  spoint_per_surfel_.resize(surfel_planes_.size());

  std::cout << "Plane type  :" << counter.transpose()
            << "; Plane number: " << surfel_planes_.size() << std::endl;

  // 给平面上色
  surfels_map_.clear();
  {
    int idx = 0;
    for (const auto &v : surfel_planes_) {
      colorPointCloudT cloud_rgb;
      pcl::copyPointCloud(v.cloud_inlier, cloud_rgb);

      size_t colorType = (idx++) % color_list_.size(); // 点云颜色
      for (auto & p : cloud_rgb) {
        p.rgba = color_list_[colorType]; // 每个点给不同颜色
      }
      surfels_map_ += cloud_rgb;
    }
  }
}

// 输入：对应的在全局地图中的点云和原始点云
// 功能：使用scan_inM找到属于所有平面的点，条件是点在平面的bbox中，且点面距离小于一定阈值，并且将对应的原始点云中的点加入进来
void SurfelAssociation::getAssociation(const VPointCloud::Ptr& scan_inM,
                                       const TPointCloud::Ptr& scan_raw,
                                       size_t selected_num_per_ring) {
  const size_t width = scan_raw->width;
  const size_t height = scan_raw->height;

  int associatedFlag[width * height];
  for (unsigned int i = 0; i < width * height; i++) {
    associatedFlag[i] = -1; // 初始化为-1，表示没有关联上
  }

  // for循环多线程
  #pragma omp parallel for num_threads(omp_get_max_threads())
  for (int plane_id = 0; plane_id < surfel_planes_.size(); plane_id++) { // 遍历所有的候选平面
    std::vector<std::vector<int>> ring_masks;
    associateScanToSurfel(plane_id, scan_inM, associated_radius_, ring_masks); //对于当前scan，遍历每个点，如果在bbox里面，且点面距离小于radius，则将其属于哪一根线加入到mask_per_ring

    // 每根线上稀疏化选selected_num_per_ring个点，以免让优化问题太庞大
    for (int h = 0; h < height; h++) {
      if (ring_masks.at(h).size() < selected_num_per_ring*2) //当前scanID上属于这个平面的点太少，不加入
        continue;
      int step = ring_masks.at(h).size() / (selected_num_per_ring + 1); // 当前线束上选中的点/3，但是根据上面的判断，每一行至少有4个点
      step = std::max(step, 1); // 选最大的，跳着选selected_num_per_ring个点，使点更均匀
      for (int selected = 0; selected < selected_num_per_ring; selected++) { // 每个ring上只选择selected_num_per_ring个点
        int w = ring_masks.at(h).at(step * (selected+1) - 1);
        associatedFlag[h*width + w] = plane_id; // 当前点对应的哪一个平面
      }
    }
  }
  // 在scanM中找到属于每个plane_id的点，每根线上只需要两个点，然后将scanM选中的点在scan_raw中找到对应的点，加入到spoint中，最后
  // in chronological order
  SurfelPoint spoint;
  for (int w = 0; w < width; w++) {
    for (int h = 0; h < height; h++) {
      if (associatedFlag[ h * width + w] == -1
          || 0 == scan_raw->at(w,h).timestamp)
        continue;
      spoint.point = Eigen::Vector3d(scan_raw->at(w,h).x,
                                     scan_raw->at(w,h).y,
                                     scan_raw->at(w,h).z);
      spoint.point_in_map = Eigen::Vector3d(scan_inM->at(w,h).x,
                                            scan_inM->at(w,h).y,
                                            scan_inM->at(w,h).z);
      spoint.plane_id = associatedFlag[h*width+w];
      spoint.timestamp = scan_raw->at(w,h).timestamp;

      spoint_per_surfel_.at(spoint.plane_id).push_back(spoint); // 每个平面对应的点云
      spoints_all_.push_back(spoint);
    }
  }
}


void SurfelAssociation::randomDownSample(int num_points_max) {
  for (const auto & v : spoint_per_surfel_) {
    if (v.size() < 20)
      continue;
    for(int i = 0; i < num_points_max; i++) {
      int random_index = rand() / (RAND_MAX) * v.size();
      spoint_downsampled_.push_back(v.at(random_index));
    }
  }
}


void SurfelAssociation::averageDownSmaple(int num_points_max) {
  for (const auto & v : spoint_per_surfel_) {
    if (v.size() < 20)
      continue;
    int d_step = v.size() / num_points_max;
    int step = d_step > 1 ? d_step : 1;
    for (int i = 0; i < v.size(); i+=step) {
      spoint_downsampled_.push_back(v.at(i));
    }
  }
}
// 默认10个点取一个，返回平面，将spoints_all_降采样，这些点时要加入优化的
void SurfelAssociation::averageTimeDownSmaple(int step) {
  for (size_t idx = 0; idx < spoints_all_.size(); idx+= step) {
    spoint_downsampled_.push_back(spoints_all_.at(idx));
  }
}
// 返回与坐标轴平行的那个轴序号
int SurfelAssociation::checkPlaneType(const Eigen::Vector3d& eigen_value,
                                      const Eigen::Matrix3d& eigen_vector,
                                      const double& p_lambda) {
  Eigen::Vector3d sorted_vec;
  Eigen::Vector3i ind;
  Eigen::sort_vec(eigen_value, sorted_vec, ind); //sorted_vec由大到小排排列，对应的序号保存在ind里面

  double p = 2*(sorted_vec[1] - sorted_vec[2]) /
             (sorted_vec[2] + sorted_vec[1] + sorted_vec[0]);

  if (p < p_lambda) { // 小于阈值，不是平面
    return -1;
  }

  int min_idx = ind[2]; // 最小特征值对应列的特征向量为平面的法向量
  Eigen::Vector3d plane_normal = eigen_vector.block<3,1>(0, min_idx); // 对应列
  plane_normal = plane_normal.array().abs(); // 先取绝对值

  Eigen::sort_vec(plane_normal, sorted_vec, ind);// 平面法向量由大到小排序，返回最小的法向量序号，返回与轴平行的序号，012代表xyz轴 //?
  return ind[2];
}

// cloud RANSAC拟合平面，平面参数保存到coeffs中，平面内点保存到cloud_inliers中，inlier点数大与20
bool SurfelAssociation::fitPlane(const VPointCloud::Ptr& cloud,
                                 Eigen::Vector4d &coeffs,
                                 VPointCloud::Ptr cloud_inliers) {
  pcl::ModelCoefficients::Ptr coefficients (new pcl::ModelCoefficients);
  pcl::PointIndices::Ptr inliers (new pcl::PointIndices);
  pcl::SACSegmentation<VPoint> seg;    /// Create the segmentation object
  // Optional
  seg.setOptimizeCoefficients (true);
  // Mandatory
  seg.setModelType (pcl::SACMODEL_PLANE);
  seg.setMethodType (pcl::SAC_RANSAC);
  seg.setDistanceThreshold (0.05);

  seg.setInputCloud (cloud);
  seg.segment (*inliers, *coefficients);

  if (inliers->indices.size () < 20) {
    return false;
  }

  for(int i = 0; i < 4; i++) {
    coeffs(i) = coefficients->values[i];
  }

  pcl::copyPointCloud<VPoint> (*cloud, *inliers, *cloud_inliers);
  return true;
}

double SurfelAssociation::point2PlaneDistance(Eigen::Vector3d &pt,
                                              Eigen::Vector4d &plane_coeff) {
  Eigen::Vector3d normal = plane_coeff.head<3>();
  double dist = pt.dot(normal) + plane_coeff(3);
  dist = dist > 0 ? dist : -dist;

  return dist;
}

// 对于当前scan，遍历每个点，如果在bbox里面，且点面距离小于radius，则将其属于哪一根线加入到mask_per_ring
void SurfelAssociation::associateScanToSurfel(
        const size_t& surfel_idx,
        const VPointCloud::Ptr& scan, // 地图中的点云
        const double& radius,
        std::vector<std::vector<int>> &ring_masks) const {
  // 获取对应平面信息
  Eigen::Vector3d box_min = surfel_planes_.at(surfel_idx).boxMin;
  Eigen::Vector3d box_max = surfel_planes_.at(surfel_idx).boxMax;
  Eigen::Vector4d plane_coeffs = surfel_planes_.at(surfel_idx).p4;

  for (int j = 0 ; j < scan->height; j++) {
    std::vector<int> mask_per_ring;
    for (int i=0; i< scan->width; i++) { // 每一行，这个点在这个平面的bbox里面，且到平面的距离小于0.05
      if (!pcl_isnan(scan->at(i,j).x) && // 点不是nan点
          scan->at(i,j).x > box_min[0] && scan->at(i,j).x < box_max[0] &&
          scan->at(i,j).y > box_min[1] && scan->at(i,j).y < box_max[1] &&
          scan->at(i,j).z > box_min[2] && scan->at(i,j).z < box_max[2]) { // 在bbox里面

          Eigen::Vector3d point(scan->at(i,j).x, scan->at(i,j).y, scan->at(i,j).z);
          if (point2PlaneDistance(point, plane_coeffs) <= radius) { // 放弃点面距离较大的点
            mask_per_ring.push_back(i);
          }
      }
    } // end of one colmun (ring)
    ring_masks.push_back(mask_per_ring); // 每一行代表这行属于当前平面的点
  }
}

}
