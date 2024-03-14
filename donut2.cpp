#include<vector>
#include<cmath>
#include<chrono>
#include<iostream>
#include<string>
#include<limits>
#include<thread>
#include<stdlib.h>


                                                                using namespace std;struct coord3D
                                                    {double x;double y;double z;coord3D(double x,double y,
                                              double z) : x(x), y(y), z(z){}bool operator==(const coord3D other
                                        )const{const double tolerance=1e-2;return abs(x-other.x)< tolerance&&abs(
                                      y-other.y) < tolerance&&abs(z - other.z)<tolerance;}coord3D operator-(const 
                                   coord3D other)const{return{x-other.x,y-other.y,z-other.z};}coord3D operator+(const 
                                 coord3D other)const{return{x+other.x,y+other.y,z+other.z};}double operator*(const coord3D
                              other)const{return x*other.x+y*other.y + z* other.z;}void normalize(){double magnitude=sqrt(x
                           *x+y*y+z*z);x/=magnitude;y /= magnitude;z /= magnitude;}coord3D rotate(const double(&matrix)[3][3]
                         )const{return{x*matrix[0][0]+y*matrix[1][0]+z*matrix[2][0],x*matrix[0][1]+y*matrix[1][1]+z*matrix[2][1]
                       ,x*matrix[0][2]+y*matrix[1][2]+z*matrix[2][2]};}};const int FPS = 20;const double R=4.0;const double r=2.2;
                      const int WINDOW_WIDTH_CHAR=61;const int WINDOW_HEIGHT_CHAR= 31;const coord3D light={0, 5, -15};const coord3D
                    view={0,0,-10};const int RAYS_PER_UNIT_X=4;const int RAYS_PER_UNIT_Y=2;const double WINDOW_HEIGHT=2.5*(R+r);const
                  double WINDOW_WIDTH=2.5*(R+r);bool point_inside_donut(coord3D p){static double _4R2=4*R*R;static double R2_minus_r2=
                 R*R-r*r;double x2_plus_y2=p.x*p.x+p.y*p.y;double inside_sq=(x2_plus_y2+p.z*p.z+R2_minus_r2);return inside_sq*inside_sq
               -_4R2*x2_plus_y2<= 0.0;}coord3D vector_intersects_donut(coord3D& loc,coord3D& dir) { static const double jump = 0.1;static
              const double jump_jump=0.01;static const double R_plus_r=R+r;static const double _R_minus_r=-R_plus_r;static const double _r
             =-r;double t1=numeric_limits<double>::infinity(),t2=numeric_limits<double>::infinity();static auto min2=[](double t,double& t1,
            double& t2){if(t<t1){t2=t1;t1=t;}else if(t<t2){t2=t;}};double temp;double x,y,z;if(dir.x!=0){temp=(R_plus_r-loc.x)/dir.x;y=loc.y
           +dir.y*temp,z=loc.z+dir.z*temp;if(y<=R_plus_r&&y>=_R_minus_r&&z<=r&&z>=_r)min2(temp,t1,t2);temp=(_R_minus_r- loc.x)/dir.x ;y=loc.y+
          dir.y*temp,z=loc.z+dir.z*temp;if(y<=R_plus_r&&y>=_R_minus_r          &&z<=r&&z>=_r)min2(temp,t1,t2);}if(dir.y!=0){temp=(R_plus_r-loc
         .y)/dir.y;x=loc.x+dir.x*temp,z=loc.z+dir.z*temp;if(x<=                     R_plus_r&&x>=_R_minus_r&&z<=r&&z>=_r)min2(temp,t1,t2);temp=
        (_R_minus_r-loc.y)/dir.y;x=loc.x+dir.x*temp,z=loc.z+                           dir.z*temp;if(x<=R_plus_r&&x>=_R_minus_r &&z<=r&& z >=_r)
       min2(temp,t1,t2);}if(dir.z!=0){temp=(r-loc.z)/dir.                                 z;x=loc.x+dir.x*temp,y=loc.y+dir.y*temp;if(x<=R_plus_r
       &&x>=_R_minus_r&&y<=R_plus_r&&y>=_R_minus_r)min2(                                   temp,t1,t2);temp=(_r-loc.z)/dir.z;x=loc.x+dir.x*temp,y
      =loc.y+dir.y*temp;if(x<=R_plus_r&&x>=_R_minus_r&&                                      y<=R_plus_r&&y>=_R_minus_r)min2(temp,t1,t2);}double t
     =t1;while(t<t2 &&!point_inside_donut({loc.x+dir.x                                        *t,loc.y+dir.y*t,loc.z+dir.z*t})){t+=jump;}if(t>=t2)
     {return{numeric_limits<double>::infinity(),/****/                                        numeric_limits<double>::infinity(),/**/numeric_limits
     <double>::infinity()};}while(point_inside_donut(                                          {loc.x+dir.x* t,loc.y+dir.y* t,loc.z+dir.z*t})){t -=
     jump_jump;}return{loc.x+dir.x*t,loc.y+dir.y *t,                                           loc.z+dir.z*t};}coord3D get_normal(coord3D& point){
     static double R2=R*R;static double _2PI=2*M_PI;                                           static double _3_over_2_PI=3*M_PI_2;double phi=atan2
     (point.y,point.x);double div = point.z/r;if(div>                                         1.0)div=1.0;else if(div<-1.0)div=-1.0;double theta =
     asin(div);bool upper=point.z>0;bool right=point.                                         x*point.x+point.y*point.y>R2;if(point.z==0.0&&right)
     {theta=0;}else if(point.z==r){theta=M_PI_2;}else                                        if(point.z==0.0&&!right){theta=M_PI;}else if(point.z==
    -r){theta=_3_over_2_PI;}else if(upper&&right){theta                                     =theta;}else if(upper&&!right){theta=M_PI-theta;}else 
    if(!upper&&!right){theta=M_PI-theta;}else{theta=_2PI+                                  theta;}double tx=-sin(phi);double ty=cos(phi);double tz
    =0;double sx=cos(phi)*-sin(theta);double sy=sin(phi)*-                               sin(theta);double sz=cos(theta);double nx=ty*sz-tz * sy;
     double ny=tz*sx-tx*sz;double nz=tx*sy-ty*sx;coord3D normal                        ={nx,ny,nz};normal.normalize();/*HI*/return normal;}float
     angle_between_vectors(coord3D& a,coord3D& b){double dot=-a.x*                 b.x-a.y*b.y-a.z*b.z;return acos(dot);}char light_to_char(double
     angle_degrees){static const string light_chars=".,-~:;=!*#$@";static    const int num_light_chars = light_chars.size();static const double
     scale=90.0/11.9;if(angle_degrees >= 90.0) return '.';int index=int((90.0-angle_degrees)/scale);return light_chars[index];/*HELLO*.*/}char*
     render_ascii_donut(double A,double B,double C){static const double x_min=-WINDOW_WIDTH_CHAR/(2*RAYS_PER_UNIT_X);static const double x_max
     =-x_min;static const double y_min=-WINDOW_HEIGHT_CHAR/(2*RAYS_PER_UNIT_Y);static const double y_max=-y_min;static const double x_inc=1.0/
     RAYS_PER_UNIT_X;static const double y_inc=1.0/RAYS_PER_UNIT_Y;static const double to_degrees=180.0/M_PI;static const int total_chars =
      WINDOW_HEIGHT_CHAR*(WINDOW_WIDTH_CHAR+1)+1;const double cosA=cos(A);const double sinA=sin(A);const double cosB=cos(B);const double sinB
       =sin(B);const double cosC=cos(C);const double sinC=sin(C);double rotation_matrix[3][3]={{cosB*cosC,sinA*sinB*cosC-cosA*sinC,cosA*sinB
       *cosC+sinA*sinC},{cosB*sinC,sinA*sinB*sinC+cosA*cosC,cosA*sinB*sinC-sinA*cosC},{-sinB,sinA*cosB,cosA*cosB}};coord3D lr=light.rotate
        (rotation_matrix);coord3D view_rotated=view.rotate(rotation_matrix);coord3D vdr=coord3D(0,0,0)-view_rotated;vdr.normalize();char*
         buffer=new char[total_chars];int i=0;for(double y=y_max;y>=y_min;y -=y_inc){for(double x =x_min;x<=x_max;x+=x_inc){double 
           inside_paren=y*cosA-view.z * /*HAYDENDIPPL*/sinA;coord3D vlr = coord3D(x, y , view.z).rotate(rotation_matrix);coord3D vip=
            vector_intersects_donut(vlr,vdr);bool inter_donut=vip.z!=numeric_limits<double>::infinity();if(!inter_donut){buffer[i++]=
            ' ';continue;}coord3D ldtd=vip-lr;ldtd.normalize();coord3D lip=vector_intersects_donut(lr,ldtd);if(lip.x==numeric_limits
              <double>::infinity()||!(vip==lip)){buffer[i++]='.';continue;}coord3D normal=get_normal(vip);double angle_degrees=
                 angle_between_vectors(normal,ldtd) * to_degrees;buffer[i++] = light_to_char(angle_degrees);}buffer[i++]='\n';}
                 buffer[i]='\0';return buffer;}chrono::milliseconds get_time(){return chrono::duration_cast<chrono::milliseconds
                 >(chrono::system_clock::now().time_since_epoch());}int main()/*DONUT*/{static const auto ms_in_frame=chrono::
                   milliseconds(1000/FPS);static const double inc_yaw=6.0*M_PI/180.0;static const double inc_pitch=2.5*M_PI
                     /180.0;static const double inc_roll=0.0*M_PI/180.0;static const double _2PI=2.0*M_PI;double yaw=0.0;
                       double pitch=0.0;double roll=0.0;while(true){chrono::milliseconds start_frame_ms = get_time();
                         char* donut=render_ascii_donut(yaw,pitch,roll);chrono::milliseconds end_render_ms=get_time
                            ();chrono::milliseconds render_duration_ms/*CODE*/=end_render_ms-start_frame_ms;chrono::
                              milliseconds sleep_duration_ms=ms_in_frame-render_duration_ms;if (sleep_duration_ms
                                  .count()>0)this_thread::sleep_for(sleep_duration_ms);system("clear");cout
                                    <<donut;yaw+=inc_yaw;pitch+=inc_pitch;roll+=inc_roll;if(yaw>_2PI)
                                          yaw-=_2PI;if(pitch>_2PI)pitch-=_2PI;if(roll>_2PI)roll-=
                                                                _2PI;}return 0;}