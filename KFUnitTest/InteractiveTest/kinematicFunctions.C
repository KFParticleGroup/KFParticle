#include <cmath>

float PT(float px,float py){
    return sqrt(px*px+py*py);
}

float P(float px,float py,float pz) {
    return sqrt(px*px+py*py+pz*pz);
}

float E(float px,float py,float pz, float mass) {
    return sqrt(px*px+py*py+pz*pz+mass*mass);
}

float Eta(float px,float py,float pz) {
    float p = P(px, py, pz);
    return 0.5*log((p+pz)/(p-pz));
}

float Phi(float px, float py){
    return atan2(py, px);
}

float Theta(float px,float py,float pz){
    return atan2(PT(px,py), pz);
}

