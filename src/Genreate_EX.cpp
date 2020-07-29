#include<iostream>
#include<fstream>
#include<iomanip>
#include<dirent.h>
#include<string>
#include<vector>
#include <iterator>
#include <regex>

using namespace std;

#define PI 3.1415926

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

/* 
   用delim指定的正则表达式将字符串in分割，返回分割后的字符串数组
   delim 分割字符串的正则表达式 
 */
std::vector<std::string> s_split(const std::string& in, const std::string& delim) {
    std::regex re{ delim };
    // 调用 std::vector::vector (InputIterator first, InputIterator last,const allocator_type& alloc = allocator_type())
    // 构造函数,完成字符串分割
    return std::vector<std::string> {
        std::sregex_token_iterator(in.begin(), in.end(), re, -1),
        std::sregex_token_iterator()
    };
}

int main(int argc, char** argv )
{
    if(argc != 2)
    {
        cout << "./LB_Ex Dir_path" <<endl;
        exit(0);
    }

    vector<string> strPrefix_ex;
    strPrefix_ex=getFiles(string(argv[1]));
	const int nCams = (int)strPrefix_ex.size();

    string outpath = "./ExPara.txt";
    ofstream out(outpath);

    for(int nc = 0; nc < nCams; nc++)
    {
        string name = string(argv[1]) + string("/") + strPrefix_ex[nc];
        ifstream in(name);
        if(!in.is_open())
        {
            cout << "Can't open " << name <<endl;
            exit(0);
        }

        string s;
        getline(in, s);
        std::vector<string> all = s_split(s, string(","));

        stringstream ss;
        double Rx, Ry, Rz, Tx, Ty, Tz;
        ss << all[0]; ss >> Tx; ss.clear(); cout << Tx << " "; 
        ss << all[1]; ss >> Ty; ss.clear(); cout << Ty << " ";
        ss << all[2]; ss >> Tz; ss.clear(); cout << Tz << " ";
        ss << all[3]; ss >> Rx; Rx = Rx / 180.0 * PI; ss.clear(); cout << Rx << " ";
        ss << all[4]; ss >> Ry; Ry = Ry / 180.0 * PI; ss.clear(); cout << Ry << " ";
        ss << all[5]; ss >> Rz; Rz = Rz / 180.0 * PI; ss.clear(); cout << Rz << " ";
        cout << endl;

        out << Rx << "  " << Ry << "  " << Rz << "  " << Tx << "  "<< Ty << "  " << Tz << endl;
        in.close();
    }

    out.close();

    return 0;
}