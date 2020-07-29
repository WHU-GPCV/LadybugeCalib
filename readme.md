
# ladybug 检校
LadybugCalib实现导出ladybug全景相机内参(等距投影+畸变参数)的功能。服务于pano-slam系统。
## 源代码
需先配置第三方依赖库：Eigen、ceres  

## 可执行程序

运行环境：ubuntu 18.04, ubuntu 16.04
运行./LBCalib   Dir_path_D2U   Dir_path_U2D

比如：
./LBCalib ../testdata/U2D_ALL_1616X1232/ ../testdata/D2U_ALL_1616X1232/

导出结果为 InnPara.txt， InvInnPara.txt

