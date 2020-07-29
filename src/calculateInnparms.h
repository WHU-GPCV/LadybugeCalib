#include "ladybugCalib.h"

#include <ceres/ceres.h>

using namespace ceres;

struct CostFunctor1 {
    CostFunctor1(double xx, double yy, double xxr, double yyr, double xx0rectified, double yy0rectified, double fInit)
        : x(xx), y(yy), xr(xxr), yr(yyr), x0rectified(xx0rectified), y0rectified(yy0rectified), 
        f(fInit), xequi(xr-x0rectified), yequi(yr-y0rectified) {}
    
    template<typename T>
    bool operator()(const T* const lamda, const T* const Ip0, const T* const k, const T* const p, const T* const c, T* residuals) const 
    { 
        T xt = T(x)-Ip0[0];
        T yt = T(y)-Ip0[1];
        T sqrtr=ceres::sqrt(xt*xt+yt*yt);
        T Rectifiedx=((T(f)*xt*ceres::tan(sqrtr/T(f)))/sqrtr);
        T Rectifiedy=((T(f)*yt*ceres::tan(sqrtr/T(f)))/sqrtr);
        T x_=lamda[0]*Rectifiedx;
        T y_=lamda[0]*Rectifiedy;
        T r2 = x_*x_ + y_*y_;
        T dx = x_*(k[0]*r2 + k[1]*r2*r2+k[2]*r2*r2*r2 + k[3]*r2*r2*r2*r2) + 2.0*p[0]*x_*y_+ p[1]*(r2+2.0*x_*x_)+ c[0]*x_ + c[1]*y_ ;
        T dy = y_*(k[0]*r2 + k[1]*r2*r2+k[2]*r2*r2*r2 + k[3]*r2*r2*r2*r2) + p[0]*(r2+2.0*y_*y_)+2.0*p[1]*x_*y_  +c[0]*y_ + c[1]*x_;
        //double xequi=xr-x0rectified;
        //double yequi=yr-y0rectified;

        residuals[0] = ceres::abs(x_-dx-T(xequi));
        residuals[1] = ceres::abs(y_-dy-T(yequi));

		return true;
    }
    
    private:
    // 观测值
    const double x;
    const double y;
    const double xr;
    const double yr;
    const double x0rectified;
    const double y0rectified;
    const double xequi;
    const double yequi;
    const double f;
};

struct CostFunctor2 {
    CostFunctor2(double xx, double yy, double xxr, double yyr, double xx0rectified, double yy0rectified, double fInit)
        : x(xx), y(yy), xr(xxr), yr(yyr), x0rectified(xx0rectified), y0rectified(yy0rectified), 
        f(fInit) {}
    
    template<typename T>
    bool operator()(const T* const lamda, const T* const Ip0, const T* const k, const T* const p, const T* const c, T* residuals) const 
    { 

        T x_ = T(xr)-T(x0rectified);
        T y_ = T(yr)-T(y0rectified);
        T r2 = x_*x_ + y_*y_;
        T dx = x_*(k[0]*r2 + k[1]*r2*r2+k[2]*r2*r2*r2 + k[3]*r2*r2*r2*r2) + 2.0*p[0]*x_*y_+ p[1]*(r2+2.0*x_*x_)+ c[0]*x_ + c[1]*y_ ;
        T dy = y_*(k[0]*r2 + k[1]*r2*r2+k[2]*r2*r2*r2 + k[3]*r2*r2*r2*r2) + p[0]*(r2+2.0*y_*y_)+2.0*p[1]*x_*y_  +c[0]*y_ + c[1]*x_;
        T xt=(x_+dx);
        T yt=(y_+dy);
        T sqrtr=ceres::sqrt(xt*xt+yt*yt);
        T Distortedx=((T(f)*xt*ceres::atan(sqrtr/(T(f))))/sqrtr);
        T Distortedy=((T(f)*yt*ceres::atan(sqrtr/(T(f))))/sqrtr);
        T xd=lamda[0]*Distortedx;
        T yd=lamda[0]*Distortedy;

        T xequi=T(x)-Ip0[0];
        T yequi=T(y)-Ip0[1];


        residuals[0]=ceres::abs(xd-xequi);
        residuals[1]=ceres::abs(yd-yequi);

		return true;
    }
    
    private:
    // 观测值
    const double x;
    const double y;
    const double xr;
    const double yr;
    const double x0rectified;
    const double y0rectified;
    const double f;
};

std::vector<std::vector<double>>  CalculateDistortPara(std::vector<std::vector<double>> points, std::vector<double> InnInitial, int nWidth_, int nHight_, std::vector<double> &Inn)
{
    int nWidth = nWidth_, nHight = nHight_; 
    double lamdaInitial, FocalLengthInitial;
    double frectified = InnInitial[4];
    double x0distorted=InnInitial[0];
    double y0distorted=InnInitial[1];
    double x0rectified=InnInitial[2];
    double y0rectified=InnInitial[3];
    std::vector<std::vector<double>> pointsForInitial;

    EstimateF(points,InnInitial,nWidth,nHight,pointsForInitial,FocalLengthInitial,lamdaInitial);

    Eigen::Matrix<double, Eigen::Dynamic, 8> Air;
    Air.resize(pointsForInitial.size()*2,8);
    Eigen::Matrix<double, Eigen::Dynamic, 1> Dir;
    Dir.resize(pointsForInitial.size()*2,1);
    //计算畸变参数初始值
    for(size_t i = 0; i<pointsForInitial.size(); i++)
    {
        double xi=pointsForInitial[i][0];
        double yi=pointsForInitial[i][1];
        double xri=pointsForInitial[i][2];
        double yri=pointsForInitial[i][3];
        double ri2=xi*xi+yi*yi;
        double ri4=ri2*ri2;
        double ri6=ri2*ri2*ri2;
        double ri8=ri2*ri2*ri2*ri2;
        Eigen::Matrix<double,1,8> Air1,Air2;
        Air1<<xi*ri2,xi*ri4,xi*ri6,xi*ri8,2*xi*yi,(ri2+2*xi*xi),xi,yi;
        Air2<<yi*ri2,yi*ri4,yi*ri6,yi*ri8,(ri2+2*yi*yi),2*xi*yi,yi,xi;
        Air.row(2*i)=Air1;
        Air.row(2*i+1)=Air2;
        Dir(2*i,0)=(xi-xri);
        Dir(2*i+1,0)=(yi-yri);
    }

    Eigen::Matrix<double, 8, 1> DistInitial=(Air.transpose()*Air).inverse()*(Air.transpose()*Dir);

    //赋初始值
    double lamda  = {lamdaInitial};
    double k[4] = {DistInitial(0), DistInitial(1), DistInitial(2), DistInitial(3)};
    double p[2] = {DistInitial(4), DistInitial(5)};
    double c[2] = {DistInitial(6), DistInitial(7)};
    double Ip0[2] = {x0rectified, y0rectified};

    Problem problem;
    for (size_t i = 0; i < points.size();++i)
    {
        problem.AddResidualBlock(
            new AutoDiffCostFunction<CostFunctor1,2,1,2,4,2,2>(
                new CostFunctor1(points[i][0],points[i][1], points[i][2],points[i][3], x0rectified, y0rectified, FocalLengthInitial)), NULL,
                &lamda,Ip0,k,p,c);
    }

    ceres::Solver::Options m_options;
    ceres::Solver::Summary m_summary;
    m_options.max_num_iterations = 25;
    m_options.linear_solver_type = ceres::DENSE_QR;
    m_options.minimizer_progress_to_stdout = true;

    Solve(m_options, &problem,&m_summary);

    Inn.clear();
    Eigen::Matrix<double, 15, 1> result;
    result<<k[0],k[1],k[2],k[3],p[0],p[1],c[0],c[1],lamda,FocalLengthInitial,Ip0[0],Ip0[1],frectified,x0rectified,y0rectified;
    for(size_t i = 0;i<15;i++)
        Inn.push_back(result(i,0));
    
	return pointsForInitial;
}

void CalculateDistortParaInv(std::vector<std::vector<double>> points, std::vector<double> InnInitial, int nWidth_, int nHight_, std::vector<double> &InnInv, std::vector<std::vector<double>> pointsForInitial)
{
    int nWidth = nWidth_, nHight = nHight_; 
    double lamdaInitial, FocalLengthInitial;
    double frectified = InnInitial[4];
    double x0distorted=InnInitial[0];
    double y0distorted=InnInitial[1];
    double x0rectified=InnInitial[2];
    double y0rectified=InnInitial[3];
    std::vector<std::vector<double>> pointsForInitial1;

    EstimateFInv(points,InnInitial,nWidth,nHight,pointsForInitial1,FocalLengthInitial,lamdaInitial);

    Eigen::Matrix<double, Eigen::Dynamic, 8> Air;
    Air.resize(pointsForInitial.size()*2,8);
    Eigen::Matrix<double, Eigen::Dynamic, 1> Dir;
    Dir.resize(pointsForInitial.size()*2,1);
    //计算畸变参数初始值
    for(size_t i = 0; i<pointsForInitial.size(); i++)
    {
        double xi=pointsForInitial[i][0];
        double yi=pointsForInitial[i][1];
        double xri=pointsForInitial[i][2];
        double yri=pointsForInitial[i][3];
        double ri2=xri*xri+yri*yri;
        double ri4=ri2*ri2;
        double ri6=ri2*ri2*ri2;
        double ri8=ri2*ri2*ri2*ri2;
        Eigen::Matrix<double,1,8> Air1,Air2;
        Air1<<xri*ri2,xri*ri4,xri*ri6,xri*ri8,2*xri*yri,(ri2+2*xri*xri),xri,yri;
        Air2<<yri*ri2,yri*ri4,yri*ri6,yri*ri8,(ri2+2*yri*yri),2*xri*yri,yri,xri;
        Air.row(2*i)=Air1;
        Air.row(2*i+1)=Air2;
        Dir(2*i,0)=(xi-xri);
        Dir(2*i+1,0)=(yi-yri);
    }

    Eigen::Matrix<double, 8, 1> DistInitial=(Air.transpose()*Air).inverse()*(Air.transpose()*Dir);

    //赋初始值
    double lamda  = {lamdaInitial};
    double k[4] = {DistInitial(0), DistInitial(1), DistInitial(2), DistInitial(3)};
    double p[2] = {DistInitial(4), DistInitial(5)};
    double c[2] = {DistInitial(6), DistInitial(7)};
    double Ip0[2] = {x0rectified, y0rectified};

    Problem problem;
    for (size_t i = 0; i < points.size();++i)
    {
        problem.AddResidualBlock(
            new AutoDiffCostFunction<CostFunctor2,2,1,2,4,2,2>(
                new CostFunctor2(points[i][0],points[i][1], points[i][2],points[i][3], x0rectified, y0rectified, FocalLengthInitial)), NULL,
                &lamda,Ip0,k,p,c);
    }

    ceres::Solver::Options m_options;
    ceres::Solver::Summary m_summary;
    m_options.max_num_iterations = 25;
    m_options.linear_solver_type = ceres::DENSE_QR;
    m_options.minimizer_progress_to_stdout = true;

    Solve(m_options, &problem,&m_summary);

    InnInv.clear();
    Eigen::Matrix<double, 15, 1> result;
    result<<k[0],k[1],k[2],k[3],p[0],p[1],c[0],c[1],lamda,FocalLengthInitial,Ip0[0],Ip0[1],frectified,x0rectified,y0rectified;
    for(size_t i = 0;i<15;i++){
        InnInv.push_back(result(i,0));
    }
}
