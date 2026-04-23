#include <random>
#include <chrono> 

#include <ctime>    // For time()
#include <cstdlib>  // For srand() and rand()

#include "AliExternalTrackParam.h"

void CreateCovMat(double cv[21]){
    srand(time(0));
    std::uniform_real_distribution<double> unif(0,20);
    std::default_random_engine re;
    re.seed(std::chrono::system_clock::now().time_since_epoch().count());

    cv[0] = unif(re);
    cv[2] = unif(re);
    cv[5] = unif(re);
    cv[9] = unif(re);
    cv[14] = unif(re);
    cv[20] = unif(re);

    cv[1] = sqrt(cv[0]*cv[2])*unif(re)/200;
    cv[3] = sqrt(cv[0]*cv[5])*unif(re)/200;
    cv[4] = sqrt(cv[2]*cv[5])*unif(re)/200;
    cv[6] = sqrt(cv[0]*cv[9])*unif(re)/200;
    cv[7] = sqrt(cv[2]*cv[9])*unif(re)/200;
    cv[8] = sqrt(cv[5]*cv[9])*unif(re)/200;
    cv[10] = sqrt(cv[0]*cv[14])*unif(re)/200;
    cv[11] = sqrt(cv[2]*cv[14])*unif(re)/200;
    cv[12] = sqrt(cv[5]*cv[14])*unif(re)/200;
    cv[13] = sqrt(cv[9]*cv[14])*unif(re)/200;
    cv[15] = sqrt(cv[0]*cv[20])*unif(re)/200;
    cv[16] = sqrt(cv[2]*cv[20])*unif(re)/200;
    cv[17] = sqrt(cv[5]*cv[20])*unif(re)/200;
    cv[18] = sqrt(cv[9]*cv[20])*unif(re)/200;
    cv[19] = sqrt(cv[14]*cv[20])*unif(re)/200;

    /*short sign = 0;
    double xyz[3] = {0,0,0}, pxpypz[3] = {0,0,0};

    AliExternalTrackParam p(xyz, pxpypz, cv, sign);

    p.GetCovarianceXYZPxPyPz(cv);*/
    
}

void PrintGlobalParams(double xyz[3], double pxpypz[3], double cv[21], short sign){
    std::cout << "parameters:\t" << xyz[0] << "\t" << xyz[1] << "\t" << xyz[2] << "\t" << pxpypz[0] << "\t" << pxpypz[1] << "\t" << pxpypz[2] << std::endl;
    std::cout << "sign:\t" << sign << std::endl;
    std::cout << "covmat:\t" << "\n" << 
    cv[0]<<"\n"<<
    cv[1]<<"\t"<<cv[2]<<"\n"<<
    cv[3]<<"\t"<<cv[4]<<"\t"<<cv[5]<<"\n"<<
    cv[6]<< "\t"<<cv[7]<<"\t"<< cv[8]<<"\t"<<cv[9]<<"\n"<<
    cv[10]<< "\t"<<cv[11]<<"\t"<<cv[12]<< "\t"<<cv[13]<<"\t"<<cv[14]<<"\n"<<
    cv[15]<< "\t"<<cv[16]<<"\t"<<cv[17]<<"\t"<<cv[18]<<"\t"<<cv[19]<<"\t"<<cv[20]<<"\n";
    cout << endl;
}

void PrintLocalParams(double x, double alpha, double params[5], double cv[15], short sign){
    std::cout << "x: " << x << "\t" << "alpha: " << alpha << std::endl;
    std::cout << "parameters:\t" << params[0] << "\t" << params[1] << "\t" << params[2] << "\t" << params[3] << "\t" << params[4] << std::endl;
    std::cout << "sign:\t" << sign << std::endl;
    std::cout << "covmat:\t" << "\n" << 
    cv[0]<<"\n"<<
    cv[1]<<"\t"<<cv[2]<<"\n"<<
    cv[3]<<"\t"<<cv[4]<<"\t"<<cv[5]<<"\n"<<
    cv[6]<< "\t"<<cv[7]<<"\t"<< cv[8]<<"\t"<<cv[9]<<"\n"<<
    cv[10]<< "\t"<<cv[11]<<"\t"<<cv[12]<< "\t"<<cv[13]<<"\t"<<cv[14]<<"\n";
    cout << endl;
}



void track_param_conversion(){

    double xyz[3] = {0.1, 0.2, 0.3};
    double pxpypz[3] = {1,4,3};

    short sign = 1;

    double cv[21] = {
        31.0647,
        -7.76618,	1.94154,
        0.107363,	-0.0268407,	1.04271,
        0.215503,	-0.0538757,	0.148689,	18.0375,
        1.43335,	-0.358338,	0.00860631,	1.71594,	14.1278,
        0.205969,	-0.0514922,	0.0170188,	0.391058,	0.179991,	2.67912
    };
    //CreateCovMat(cv);

    cout << "\n -=-=-=-  Generated Global Params -=-=-=-=-\n";
    PrintGlobalParams(xyz,pxpypz, cv, sign);
    cout << endl;

    AliExternalTrackParam p(xyz, pxpypz, cv, sign);

    double glob_xyz[3], glob_pxpypz[3];
    double glob_covmat[21];

    p.GetXYZ(glob_xyz);
    p.GetPxPyPz(glob_pxpypz);
    cout << p.GetCovarianceXYZPxPyPz(glob_covmat) << endl;

    cout << "\n -=-=-=-  AliExt.. global Params - extracted -=-=-=-=-\n";
    PrintGlobalParams(glob_xyz, glob_pxpypz, glob_covmat, p.Charge());
    

}