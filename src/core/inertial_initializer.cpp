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
#include <core/inertial_initializer.h>
#include <utils/math_utils.h>

namespace licalib {

bool InertialInitializer::EstimateRotation(
        TrajectoryManager::Ptr traj_manager,
        const Eigen::aligned_vector<LiDAROdometry::OdomData>& odom_data) {

  int flags = kontiki::trajectories::EvalOrientation;
  std::shared_ptr<kontiki::trajectories::SplitTrajectory> p_traj
          = traj_manager->getTrajectory(); // kontiki优化结果中获取轨迹

  Eigen::aligned_vector<Eigen::Matrix4d> A_vec;
  // 每次取出相邻两个数据
  for (size_t j = 1; j < odom_data.size(); ++j) {
    size_t i = j - 1;
    double ti = odom_data.at(i).timestamp;
    double tj = odom_data.at(j).timestamp;
    if (tj >= p_traj->MaxTime()) // 后一个里程计数据超出时间范围，则跳出
      break;
    auto result_i = p_traj->Evaluate(ti, flags); // 获取对应时刻IMU的ti时刻的姿态
    auto result_j = p_traj->Evaluate(tj, flags);
    Eigen::Quaterniond delta_qij_imu = result_i->orientation.conjugate()
                                       * result_j->orientation;

    Eigen::Matrix3d R_Si_toS0 = odom_data.at(i).pose.topLeftCorner<3,3>();
    Eigen::Matrix3d R_Sj_toS0 = odom_data.at(j).pose.topLeftCorner<3,3>();
    Eigen::Matrix3d delta_ij_sensor = R_Si_toS0.transpose() * R_Sj_toS0; // dT
    Eigen::Quaterniond delta_qij_sensor(delta_ij_sensor);

    // 求解超参
    Eigen::AngleAxisd R_vector1(delta_qij_sensor.toRotationMatrix()); // 转换为轴角形式
    Eigen::AngleAxisd R_vector2(delta_qij_imu.toRotationMatrix());
    double delta_angle = 180 / M_PI * std::fabs(R_vector1.angle() - R_vector2.angle()); // 旋转绝对角度的差值
    double huber = delta_angle > 1.0 ? 1.0/delta_angle : 1.0; // 两个旋转角度差的阈值1.0rad

    Eigen::Matrix4d lq_mat = mathutils::LeftQuatMatrix(delta_qij_sensor); // 构建左乘和右乘矩阵
    Eigen::Matrix4d rq_mat = mathutils::RightQuatMatrix(delta_qij_imu);
    A_vec.push_back(huber * (lq_mat - rq_mat));
  }
  size_t valid_size = A_vec.size();
  if (valid_size < 15) {
    return false;
  }
  Eigen::MatrixXd A(valid_size * 4, 4);
  for (size_t i = 0; i < valid_size; ++i)
    A.block<4, 4>(i * 4, 0) = A_vec.at(i);

  Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV); // 奇异值从大到小排列
  Eigen::Matrix<double, 4, 1> x = svd.matrixV().col(3);
  Eigen::Quaterniond q_ItoS_est(x); // vector4d初始化，实部在后
  // 方式一：4个标量
  // Quaterniond q1(1, 2, 3, 4);   // 第一种方式  实部为1 ，虚部234
  // 方式二：  Vector4d
  // Quaterniond q2(Vector4d(1, 2, 3, 4));  // 第二种方式  实部为4 ，虚部123

  Eigen::Vector4d cov = svd.singularValues();
  std::cout<<" singular"<<svd.singularValues()<<std::endl;

  if (cov(2) > 0.25) {
    q_ItoS_est_ = q_ItoS_est;
    rotaion_initialized_ = true;
    return true;
  } else {
    return false;
  }
}

}
