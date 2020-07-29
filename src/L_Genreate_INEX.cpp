#include<iostream>
#include<fstream>
#include<iomanip>
#include<dirent.h>
#include<string>

//#include "calibGenerateIn.h"
#include "calculateInnparms.h"

using namespace std;

vector<string> getFiles(string cate_dir)
{
    vector<string> files;//存放文件名

    DIR *dir;
    struct dirent *ptr;

    if ((dir=opendir(cate_dir.c_str())) == NULL)
    {
        perror("Open dir error...");
        exit(1);
    }

    while ((ptr=readdir(dir)) != NULL)
    {
        if(strcmp(ptr->d_name,".")==0 || strcmp(ptr->d_name,"..")==0)    ///current dir OR parrent dir
                continue;
        else if(ptr->d_type == 8)    ///file
            files.push_back(ptr->d_name);
        else if(ptr->d_type == 10)    ///link file
            continue;
        else if(ptr->d_type == 4)    ///dir
        {
            files.push_back(ptr->d_name);
        }
    }
    closedir(dir);

    //排序，按从小到大排序
    sort(files.begin(), files.end());
    return files;
}

int main(int argc, char **argv)
{
	if(argc != 3)
	{
		cout << "./LBCalib   Dir_path_D2U   Dir_path_U2D" << endl;
		exit(0);
	}

	vector<string> strPrefix_d2u;
    strPrefix_d2u=getFiles(string(argv[1]));
	const int nCams = (int)strPrefix_d2u.size();

	vector<string> strPrefix_u2d;
    strPrefix_u2d=getFiles(string(argv[2]));

	if((int)strPrefix_u2d.size() != nCams)
	{
		cout << " the number of d2u is  not equal with the number of u2d!!!  " << endl;
		exit(0);
	}

    string Innoutpath = "./InnPara.txt";
    string InvInnoutpath = "./InvInnPara.txt";

	// optimizer
	cout << endl << "===========start optimization of Inn-InvInn Parameters =====================" << endl;
	google::InitGoogleLogging(argv[0]);

	std::vector<std::vector<double>> Inn(nCams);
	std::vector<std::vector<double>> InnInv(nCams);
	std::ofstream out1(Innoutpath);
	std::ofstream out2(InvInnoutpath);
	for (int i = 0; i < nCams; i++)
	{
		int nWidth, nHight;
		nWidth = 1616;
		nHight = 1232;

		// For D2U
		string path_D2U = string(argv[1]) + "/" + strPrefix_d2u[i];
		std::ifstream in1(path_D2U);
		if(!in1.is_open())
		{
			cout << "Can't open" << path_D2U << endl;
			exit(0);
		}

		cout << "Reading " << path_D2U <<endl;

		string s; 
		getline(in1, s);

		std::vector<std::vector<std::pair<double, double>>> maps_d2u;
		for(int y = 0; y < nHight; y++)
		{
			std::vector<std::pair<double, double>> t;
			for(int x = 0; x< nWidth; x++)
			{
				getline(in1, s);
				stringstream ss;
				ss << s;
				int v1, v2; double v3, v4;
				ss >> v1; ss >> v2;
				if(v1 != y || v2 != x){
					cout<<"LLLLLLLLLLLLLLLLLLL"<<endl;
					exit(0);
				}

				ss >> v3; ss >> v4;
				std::pair<double, double> tt(v3, v4);
				t.push_back(tt);
			}
			maps_d2u.push_back(t);
		}

		in1.close();

		// For U2D
		string path_U2D = string(argv[2]) + "/" + strPrefix_u2d[i];
		std::ifstream in2(path_U2D);
		if(!in2.is_open())
		{
			cout << "Can't open" << path_U2D << endl;
			exit(0);
		}

		cout << "Reading " << path_U2D <<endl;

		getline(in2, s);

		std::vector<std::vector<std::pair<double, double>>> maps_u2d;
		for(int y = 0; y < nHight; y++)
		{
			std::vector<std::pair<double, double>> t;
			for(int x = 0; x< nWidth; x++)
			{
				getline(in2, s);
				stringstream ss;
				ss << s;
				int v1, v2; double v3, v4;
				ss >> v1; ss >> v2;
				if(v1 != y || v2 != x){
					cout<<"LLLLLLLLLLLLLLLLLLL"<<endl;
					exit(0);
				}

				ss >> v3; ss >> v4;
				std::pair<double, double> tt(v3, v4);
				t.push_back(tt);
			}
			maps_u2d.push_back(t);
		}

		in2.close();

		// Input InitParam
		cout << "Please input x0rectified y0rectified frectified:" <<endl;
		double x0distorted, y0distorted, x0rectified, y0rectified, frectified;
		cin >> x0rectified >> y0rectified >> frectified;
		cout << endl;
		int x0r_i = (int)round(x0rectified); int y0r_i = (int)round(y0rectified); 
		std::pair<double, double> ii =  maps_u2d[y0r_i][x0r_i];
		x0distorted = ii.second; y0distorted = ii.first;

		vector<double> InnInitial;
		InnInitial.resize(6);
		InnInitial[0] = x0distorted;
		InnInitial[1] = y0distorted;
		InnInitial[2] = x0rectified;
		InnInitial[3] = y0rectified;
		InnInitial[4] = frectified;

		// Start  Optmize
		cout << "optmize " << i + 1 << " camera!" << endl;
		std::vector<std::vector<double>> points;
		for (int x = 57; x < nWidth; x = x + 150)
		{
			for (int y = 75; y < nHight; y = y + 120)
			{
				double Rectifiedx, Rectifiedy;
				std::pair<double, double> ii =  maps_d2u[y][x];
				Rectifiedx = ii.second;
				Rectifiedy = ii.first;

				std::vector<double> drp;
				drp.resize(4);
				drp[0] = x;
				drp[1] = y;
				drp[2] = Rectifiedx;
				drp[3] = Rectifiedy;
				
				points.emplace_back(drp);
			}
		}

		cout << "Inn ... ..." << endl;
		std::vector<std::vector<double>> pointsForInitial = CalculateDistortPara(points, InnInitial, nWidth, nHight, Inn[i]);
		cout << "InvInn ... ..." << endl;
		CalculateDistortParaInv(points, InnInitial, nWidth, nHight, InnInv[i], pointsForInitial);

		out1 << i << " ";
		for (size_t n = 0; n < Inn[i].size(); n++) {

			if (i != Inn[n].size() - 1) {
				out1 << Inn[i][n] << " ";
			}
			else {
				out1 << Inn[i][n] << std::endl;
			}
		}
		out1 << endl;

		out2 << i << " ";
		for (size_t n = 0; n < InnInv[i].size(); n++) {

			if (i != InnInv[i].size() - 1) {
				out2 << InnInv[i][n] << " ";
			}
			else {
				out2 << InnInv[i][n] << std::endl;
			}
		}
		out2 << endl;

	}

	out1.close();
	out2.close();
	cout <<endl << "============================================="<< endl 
		<< "Inn , InvInn and Ex parameters output!" << endl << endl;

    return 0;
}