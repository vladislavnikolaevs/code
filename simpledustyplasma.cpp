#include <iostream>
#include <cmath>
#include <stdio.h>
#include <cstdlib>
#include <ctime>

const int N = 2;
const double boltzmann = 1.38e-16;
const double m=1e-10;
const double  dt=1e-4;
const double frict = 20.0;

struct kinematic {
  double x,y,z;
 };

kinematic call_norm_distr() {
   double s,x,y;
   kinematic z;
   //If 0<s<=1, the value is generated (Box-Muller method)
   s = 1.5;
   while (s > 1.0) {
     x = 0.0;
     y = 0.0;
     while ((pow(x,2.0) + pow(y,2.0)) == 0.0) {
      x = -1.0 + rand() % 20001 / 10000.0;
      y = -1.0 + rand() % 20001 / 10000.0;
     }
     s = pow(x,2.0) + pow(y,2.0);
   }
   z.x = x*sqrt(-2*log(s)/s);
   z.y = y*sqrt(-2*log(s)/s);
   s = 1.5;
   while (s > 1.0) {
    x = 0.0;
    y = 0.0;
    while ((pow(x,2.0) + pow(y,2.0)) == 0.0) {
     x = -1.0 + rand() % 20001 / 10000.0;
     y = -1.0 + rand() % 20001 / 10000.0;
    }
    s = pow(x,2.0) + pow(y,2.0);
   }
   z.z = x*sqrt(-2*log(s)/s);
   return z;
 };

class Particles {
  public:
  int num;
  kinematic r,d,v,a,b,random;
  kinematic rn[N];
  double E,dx,dy,dz,r1,r2,c,l1,l2,q,temp;

  double distance(kinematic r, kinematic d) {
    l2 = (r.x-d.x)*(r.x-d.x)+(r.y-d.y)*(r.y-d.y)+(r.z-d.z)*(r.z-d.z);
	  l1 = sqrt(l2);
    return l1;
  }

  double energy_single (kinematic r, kinematic v, kinematic rn[],int num, double q, double alpha, double betta, double kappa) {
	   E = 0.5*m*(pow(v.x,2)+pow(v.y,2)+pow(v.z,2))+0.5*alpha*q*(pow(r.x,2)+pow(r.y,2))+0.5*betta*q*pow(r.z,2);
      for (int i = 0;i<N;i++)  {
        if(i!=num) {
    	     dx = r.x-rn[i].x;
           dy = r.y-rn[i].y;
           dz = r.z-rn[i].z;
           r2=dx*dx+dy*dy+dz*dz;
           r1=sqrt(r2);
           E += 0.5*pow(q,2.)*exp(-kappa*r1)/r1; // добавил 0.5 чтобы не учитывать энергию взаимодействия частиц два раза
    		}
    		else {
      	   E+=0.;
    		}
  	}
   	return E;
  }

  kinematic position_new (kinematic r, kinematic v, kinematic a, double dt) { //I am counting the new coordinate of a single particle
   d.x = r.x + v.x*dt + 0.5*a.x*pow(dt,2);
   d.y = r.y + v.y*dt + 0.5*a.y*pow(dt,2);
   d.z = r.z + v.z*dt + 0.5*a.z*pow(dt,2);
   return d;
  }

  kinematic acceleration_a (kinematic r, kinematic rn[], int num, double q, double alpha, double betta, double kappa, kinematic random, double temp) { //I am counting the new acceleration of a single particle
    a.x = -alpha*q*r.x/m - frict*v.x + sqrt(2.0*m*frict*boltzmann*temp/dt)*random.x/m;
    a.y = -alpha*q*r.y/m - frict*v.y + sqrt(2.0*m*frict*boltzmann*temp/dt)*random.y/m;
    a.z = -betta*q*r.z/m - frict*v.z + sqrt(2.0*m*frict*boltzmann*temp/dt)*random.z/m;
    for (int i = 0;i<N;i++) {
      if(i!=num) {
        dx = r.x-rn[i].x;
        dy = r.y-rn[i].y;
        dz = r.z-rn[i].z;
        r2=dx*dx+dy*dy+dz*dz;
        r1=sqrt(r2);
        c=pow(q,2)*exp(-kappa*r1)*(1/r1 + kappa)/r2/m;
        a.x+=dx*c;
        a.y+=dy*c;
        a.z+=dz*c;
    	}
    	else {
        a.x+=0.;
        a.y+=0.;
        a.z+=0.;
    	}
  	}
    return a;
  }

  kinematic acceleration_b (kinematic r, kinematic rn[], int num, double q, double alpha, double betta, double kappa, kinematic random, double temp) { //I am counting the new acceleration of a single particle
    b.x = -alpha*q*r.x/m - frict*v.x;// + sqrt(2.0*m*frict*boltzmann*temp/dt)*random.x/m;
    b.y = -alpha*q*r.y/m - frict*v.y;// + sqrt(2.0*m*frict*boltzmann*temp/dt)*random.y/m;
    b.z = -betta*q*r.z/m - frict*v.z;// + sqrt(2.0*m*frict*boltzmann*temp/dt)*random.z/m;
  	for (int i = 0;i<N;i++) {
      if(i!=num) {
        dx = r.x-rn[i].x;
        dy = r.y-rn[i].y;
        dz = r.z-rn[i].z;
        r2=dx*dx+dy*dy+dz*dz;
        r1=sqrt(r2);
        c=pow(q,2)*exp(-kappa*r1)*(1/r1 + kappa)/r2/m;
        b.x+=dx*c;
        b.y+=dy*c;
        b.z+=dz*c;
      }
      else {
        b.x+=0.;
        b.y+=0.;
        b.z+=0.;
      }
    }
    return b;
  }

  kinematic velocity_new (kinematic a, kinematic b, double dt) { //I am counting the new velocity of a single particle
    v.x += 0.5*(a.x+b.x)*dt;
    v.y += 0.5*(a.y+b.y)*dt;
    v.z += 0.5*(a.z+b.z)*dt;
    return v;
  }
};

int main(int argc, char *argv[]) {
  Particles dust[N];
  int i,j,k,l,temp_count,temp_num=60;
  double t,Q,temp,current,density,alpha,betta,kappa,q;
  double ambi = 4.5e-2;
  kinematic h[N];
  double distance[N][N], dist[N*N],inter_dist_temp;
  FILE *distancefile;
  srand(time(NULL));


  //I am opening the file for the pair corellation function data
  distancefile = fopen("debyepotentialdistance.txt","w");
  if (distancefile == NULL) {
    printf("File failed to open \n");
    exit(1);
  }

  
	density = 1e8;
	current = 2.27e-3;

      for (temp_count=0;temp_count<temp_num;temp_count++) {
        temp = 50.0 + 4.0*temp_count;
        //Here I am defining the particle's charge using my own calculations from Mayorov's article
 	q = 4.8e-10*(-194.859 + 1914.67*pow(temp,0.0951685));
        //Here I am defining screening length in dusty plasma
        kappa = sqrt(4*3.141592*density*4.8*4.8*1e-20/boltzmann/temp);
        //This is the trap
        alpha = ambi;
        betta = 15*alpha; //for the two-dimensional layer

        int N2 = (int) N/2.0;
        for(i=0;i<N;i++) { //I am creating a set of coordinates and velocities for particles
          if (i < N2) {
            dust[i].r.x = (rand() % 40)/1000.;
            dust[i].r.y = (rand() % 70)/1000.;
            dust[i].r.z = 0.09;
            for (j=0;j<i;j++) {
              if (dust[i].r.x == dust[j].r.x && dust[i].r.y == dust[j].r.y && dust[i].r.z == dust[j].r.z) {
                dust[i].r.x += (rand() % 30000) / 100000.;
                dust[i].r.y += (rand() % 30000) / 100000.;
                dust[i].r.z += (rand() % 30000) / 100000.;
              }
            }
          }
          else {
            dust[i].r.x = -(rand() % 40)/1000.;
            dust[i].r.y = -(rand() % 80)/1000.;
            dust[i].r.z = -0.09;
            for (j=0;j<i;j++) {
              if (dust[i].r.x == dust[j].r.x && dust[i].r.y == dust[j].r.y && dust[i].r.z == dust[j].r.z) {
                dust[i].r.x -= (rand() % 30000) / 100000.;
                dust[i].r.y -= (rand() % 30000) / 100000.;
                dust[i].r.z -= (rand() % 30000) / 100000.;
              }
            }
          }
          dust[i].v.x = 0.;
          dust[i].v.y = 0.;
          dust[i].v.z = 0.;
        }

        t = 0.;
        int z = 2e4; //duration of the run

        inter_dist_temp = 0.0;

        for(int probeg=0;probeg<200;probeg++) {

        for(i=0;i<z;i++) { //I am starting the run

          // I am counting INTER-PARTICLE DISTANCES
          if (i == z-1) {
            for (j=0;j<N;j++) {
              for (k=0;k<N;k++) {
                if (k!=j) {
                  dust[j].distance(dust[j].r,dust[k].r);
                  distance[j][k] = dust[j].l1;
                  //fprintf(correlation,"%lf",dust[j].l1);
                  }
                else {
                  distance[j][k] = 0.0;
                }
              }
            }
          }

          // I am defining the h[N] auxiliary array of particle coordinates
          for(j=0;j<N;j++) {
            h[j].x=dust[j].r.x;
            h[j].y=dust[j].r.y;
            h[j].z=dust[j].r.z;
          }

          // I am counting the total energy Q of my system (kinetical + potential only once)
          Q=0.;
          for(j=0;j<N;j++) {
            dust[j].energy_single(dust[j].r,dust[j].v,h,j,q,alpha,betta,kappa);
            Q += dust[j].E*1e9;
          }

          // I am counting accelerations for the current position
          for(j=0;j<N;j++) {
            dust[j].random = call_norm_distr();
            dust[j].acceleration_a(dust[j].r,h,j,q,alpha,betta,kappa,dust[j].random,temp);
          }

          // I am counting new position
          for(j=0;j<N;j++) {
            dust[j].position_new(dust[j].r,dust[j].v,dust[j].a,dt);
          }

          // I am redefining the h[N] auxiliary array
          for(j=0;j<N;j++) {
            h[j].x=dust[j].d.x;
            h[j].y=dust[j].d.y;
            h[j].z=dust[j].d.z;
          }

          // I am counting acceleration for the new positions
          for(j=0;j<N;j++) {
            dust[j].random = call_norm_distr();
            dust[j].acceleration_b(dust[j].d,h,j,q,alpha,betta,kappa,dust[j].random,temp);
          }
          //printf("Second check b %lf %lf %lf \n", dust[0].b.x, dust[0].b.y,dust[0].b.z);

          // I am counting new velocities
          for(j=0;j<N;j++) {
            dust[j].velocity_new(dust[j].a,dust[j].b,dt);
          }

          // Now I am redefining current positions
          for(j=0;j<N;j++) {
            dust[j].r.x=dust[j].d.x;
            dust[j].r.y=dust[j].d.y;
            dust[j].r.z=dust[j].d.z;
          }

          t+=dt;
        }
        
        // I am transforming distances 2d array into 1d array
        for (j = 0;j<N;j++) {
          for (k=0;k<N;k++) {
            dist[j*N + k]=distance[k][j];
          }
        }

        //I am arranging distances' 1d array
        double disth;
        for (j=0; j<N*N; j++) {
          for (k=j+1; k<N*N; k++) {
            if (dist[j]>dist[k]) {
              disth = dist[j];
              dist[j] = dist[k];
              dist[k] = disth;
            }
          }
        }

        inter_dist_temp += 1/200.0*(dist[0]+dist[1]+dist[2]+dist[3]);
      }
	fprintf(distancefile,"%lf %lf \n",temp,inter_dist_temp);

      }

  return 0;
  }
