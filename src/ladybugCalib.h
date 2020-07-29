#include <iostream>
#include <vector>
#include <string>
#include <math.h>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>


inline Eigen::Matrix<double, 1, 2> Rectifiedxy2Distortedxy(Eigen::Matrix<double, 1, 2> Rectifiedxy, double f)
{
    double xr=Rectifiedxy(0,0);
    double yr=Rectifiedxy(0,1);
    double sqrtr=std::sqrt(xr*xr+yr*yr);
    double Distortedx=((f*xr*std::atan(sqrtr/f))/(sqrtr));
    double Distortedy=((f*yr*std::atan(sqrtr/f))/(sqrtr));
    Eigen::Matrix<double, 1, 2> Distortedxy;
    Distortedxy(0,0) = Distortedx;
    Distortedxy(0,1) = Distortedy;

    return Distortedxy;
}

inline Eigen::Matrix<double, 1, 2> Distortedxy2Rectifiedxy(Eigen::Matrix<double, 1, 2> Distortedxy, double f)
{
    double x_=Distortedxy(0,0);
    double y_=Distortedxy(0,1);
    double sqrtr=std::sqrt(x_*x_+y_*y_);
    double Rectifiedx=((f*x_*std::tan(sqrtr/f))/sqrtr);
    double Rectifiedy=((f*y_*std::tan(sqrtr/f))/sqrtr);
    Eigen::Matrix<double, 1, 2> Rectifiedxy;
    Rectifiedxy(0,0) = Rectifiedx;
    Rectifiedxy(0,1) = Rectifiedy;

    return Rectifiedxy;
}

void EstimateF(std::vector<std::vector<double>> points, std::vector<double> InnInitial, int nWidth, int nHight, std::vector<std::vector<double>> &pointsForInitial, double &f, double &lamda)
{
    double x0distorted=InnInitial[0];
    double y0distorted=InnInitial[1];
    double x0rectified=InnInitial[2];
    double y0rectified=InnInitial[3];
    Eigen::Matrix<double, 1, 2> xy0rectified;
    xy0rectified<<x0rectified,y0rectified;

    pointsForInitial.resize(points.size());

    double minsum;
    minsum=9999999;
    // 选出合适的f值
    for(double ff = 910; ff<990; ff = ff+0.1)
    {
        Eigen::Matrix<double, 4, 2> DistortedConner;
        Eigen::Matrix<double, 4, 2> RectifiedConnerC;
        DistortedConner<<0,0,nWidth-1,0,0,nHight-1,nWidth-1,nHight-1;
        
        for(int i = 0; i<4; i++)
        {
            DistortedConner.row(i) = DistortedConner.row(i)-xy0rectified;
            RectifiedConnerC.row(i) = Distortedxy2Rectifiedxy(DistortedConner.row(i),ff);
        }
        
        double maxXDist=std::max(DistortedConner(0,0),std::max(DistortedConner(1,0),std::max(DistortedConner(2,0),DistortedConner(3,0))));
        double minXDist=std::min(DistortedConner(0,0),std::min(DistortedConner(1,0),std::min(DistortedConner(2,0),DistortedConner(3,0))));
        double maxYDist=std::max(DistortedConner(0,1),std::max(DistortedConner(1,1),std::max(DistortedConner(2,1),DistortedConner(3,1))));
        double minYDist=std::min(DistortedConner(0,1),std::min(DistortedConner(1,1),std::min(DistortedConner(2,1),DistortedConner(3,1))));

        double maxXRect=std::max(RectifiedConnerC(0,0),std::max(RectifiedConnerC(1,0),std::max(RectifiedConnerC(2,0),RectifiedConnerC(3,0))));
        double minXRect=std::min(RectifiedConnerC(0,0),std::min(RectifiedConnerC(1,0),std::min(RectifiedConnerC(2,0),RectifiedConnerC(3,0))));
        double maxYRect=std::max(RectifiedConnerC(0,1),std::max(RectifiedConnerC(1,1),std::max(RectifiedConnerC(2,1),RectifiedConnerC(3,1))));
        double minYRect=std::min(RectifiedConnerC(0,1),std::min(RectifiedConnerC(1,1),std::min(RectifiedConnerC(2,1),RectifiedConnerC(3,1))));
        int RectImgWidth=std::ceil(maxXRect-minXRect);
        int RectImgHeight=std::ceil(maxYRect-minYRect);
        double dispt=0;

		Eigen::Matrix<double, 1, 2> Distortedxy, RectifiedxyCal, Rectifiedxy;

        for(size_t i = 0;i<points.size();i++)
        {
            
            Distortedxy<<points[i][0]-x0rectified, points[i][1]-y0rectified;
            RectifiedxyCal=Distortedxy2Rectifiedxy(Distortedxy,ff);
            RectifiedxyCal=RectifiedxyCal*1.0*nWidth/RectImgWidth;
            Rectifiedxy<<points[i][2]-x0rectified,points[i][3]-y0rectified;
            dispt=dispt+(((Rectifiedxy-RectifiedxyCal).cwiseAbs()).sum());
        }
        
        if(dispt<minsum)
        {
            minsum=dispt;
            f=ff;
            lamda = 1.0*nWidth/RectImgWidth;
        }    
    }

    for(size_t i = 0;i<points.size();i++)
    {
        Eigen::Matrix<double, 1, 2> Distortedxy,RectifiedxyCal,Rectifiedxy;
        Distortedxy<<points[i][0]-x0rectified, points[i][1]-y0rectified;
        RectifiedxyCal=Distortedxy2Rectifiedxy(Distortedxy,f);
        RectifiedxyCal=RectifiedxyCal*lamda;
        Rectifiedxy<<points[i][2]-x0rectified,points[i][3]-y0rectified;
        std::vector<double> pfi_(4);
        pfi_[0] = RectifiedxyCal(0,0); pfi_[1] = RectifiedxyCal(0,1);
        pfi_[2] = Rectifiedxy(0,0); pfi_[3] = Rectifiedxy(0,1);
        pointsForInitial[i]=pfi_;
    }

}

void EstimateFInv(std::vector<std::vector<double>> points, std::vector<double> InnInitial, int nWidth, int nHight, std::vector<std::vector<double>> &pointsForInitial, double &f, double &lamda)
{
    double x0distorted=InnInitial[0];
    double y0distorted=InnInitial[1];
    double x0rectified=InnInitial[2];
    double y0rectified=InnInitial[3];
    Eigen::Matrix<double, 1, 2> xy0rectified;
    xy0rectified<<x0rectified,y0rectified;

    pointsForInitial.resize(points.size());

    //Eigen::Matrix<double, Eigen::Dynamic, 1> Distortedx=points.col(0);
    //Eigen::Matrix<double, Eigen::Dynamic, 1> Distortedy=points.col(1);
    //Eigen::Matrix<double, Eigen::Dynamic, 1> Rectifiedx=points.col(2);
    //Eigen::Matrix<double, Eigen::Dynamic, 1> Rectifiedy=points.col(3);

    double minsum;
    minsum=9999999;

    // 选出合适的f值
    for(double ff = 390; ff<410; ff = ff+0.1)
    {
        Eigen::Matrix<double, 4, 2> RectifiedConner;
        Eigen::Matrix<double, 4, 2> DistortedConnerC;
        RectifiedConner<<0,0,nWidth-1,0,0,nHight-1,nWidth-1,nHight-1;
        
        for(int i = 0; i<4; i++)
        {
            RectifiedConner.row(i) = RectifiedConner.row(i)-xy0rectified;
            DistortedConnerC.row(i) = Rectifiedxy2Distortedxy(RectifiedConner.row(i),ff);
        }
        
        double maxXDist=std::max(DistortedConnerC(0,0),std::max(DistortedConnerC(1,0),std::max(DistortedConnerC(2,0),DistortedConnerC(3,0))));
        double minXDist=std::min(DistortedConnerC(0,0),std::min(DistortedConnerC(1,0),std::min(DistortedConnerC(2,0),DistortedConnerC(3,0))));
        double maxYDist=std::max(DistortedConnerC(0,1),std::max(DistortedConnerC(1,1),std::max(DistortedConnerC(2,1),DistortedConnerC(3,1))));
        double minYDist=std::min(DistortedConnerC(0,1),std::min(DistortedConnerC(1,1),std::min(DistortedConnerC(2,1),DistortedConnerC(3,1))));

        double maxXRect=std::max(RectifiedConner(0,0),std::max(RectifiedConner(1,0),std::max(RectifiedConner(2,0),RectifiedConner(3,0))));
        double minXRect=std::min(RectifiedConner(0,0),std::min(RectifiedConner(1,0),std::min(RectifiedConner(2,0),RectifiedConner(3,0))));
        double maxYRect=std::max(RectifiedConner(0,1),std::max(RectifiedConner(1,1),std::max(RectifiedConner(2,1),RectifiedConner(3,1))));
        double minYRect=std::min(RectifiedConner(0,1),std::min(RectifiedConner(1,1),std::min(RectifiedConner(2,1),RectifiedConner(3,1))));
        int DistortedImgWidth=std::ceil(maxXDist-minXDist);
        int DistortedImgHeight=std::ceil(maxYDist-minYDist);
        double dispt=0;
        for(size_t i = 0;i<points.size();i++)
        {
            Eigen::Matrix<double, 1, 2> Rectifiedxy,DistortedxyCal,Distortedxy;
            Rectifiedxy<<points[i][2]-x0rectified, points[i][3]-y0rectified;
            DistortedxyCal=Rectifiedxy2Distortedxy(Rectifiedxy,ff);
            DistortedxyCal=DistortedxyCal*1.0*nWidth/DistortedImgWidth;
            Distortedxy<< points[i][0]-x0rectified, points[i][1]-y0rectified;
            dispt=dispt+(((Distortedxy-DistortedxyCal).cwiseAbs()).sum());
        }
        
        if(dispt<minsum)
        {
            minsum=dispt;
            f=ff;
            lamda=1.0*nWidth/DistortedImgWidth;
        }    
    }

    for(size_t i = 0;i<points.size();i++)
    {
        Eigen::Matrix<double, 1, 2> Rectifiedxy,RectifiedxyCal,Distortedxy;
        Distortedxy<<points[i][0]-x0rectified, points[i][1]-y0rectified;
        RectifiedxyCal=Distortedxy2Rectifiedxy(Distortedxy,f);
        RectifiedxyCal=RectifiedxyCal*lamda;
        Rectifiedxy<<points[i][2]-x0rectified,points[i][3]-y0rectified;
        std::vector<double> pfi_(4);
        pfi_[0] = RectifiedxyCal(0,0); pfi_[1] = RectifiedxyCal(0,1);
        pfi_[2] = Rectifiedxy(0,0); pfi_[3] = Rectifiedxy(0,1);
        pointsForInitial[i]=pfi_;
    }
}

