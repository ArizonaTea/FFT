// The class of complex numbers (complex)

#ifndef _complex
#define _complex

#include <iostream>
#include <math.h>

using namespace std;

class complex{
		double re,im;
	public:
   	complex (double=0, double=0);

      friend complex operator+ (complex, complex);
      friend complex operator- (complex,complex);
      friend complex operator* (complex,complex);
      friend complex operator/ (complex,complex);

      complex & operator+= (complex);
      complex & operator-= (complex);
      complex & operator*= (complex);
      complex & operator/= (complex);
      complex & operator= (complex);

      friend complex operator+ (complex);
      friend complex operator- (complex);

      friend int operator== (complex,complex);
      friend int operator!= (complex,complex);

      friend double real (complex);
      friend double imag (complex);
      friend complex conj (complex);
      friend double arg (complex);
      friend double norm (complex);
      friend double abs (complex);
      friend complex polar (double,double);

      friend complex exp (complex);
      friend complex log (complex);
      friend complex log10 (complex);
      friend complex sqrt (complex);
      friend complex power (const complex,int);
      friend complex power (const complex,double);
      friend complex power (double,const complex);
      friend complex power (const complex,const complex);
      friend complex sin (complex);
      friend complex cos (complex);
      friend complex tan (complex);
      friend complex sinh (complex);
      friend complex cosh (complex);
      friend complex tanh (complex);
      friend complex asin (complex);
      friend complex acos (complex);
      friend complex atan (complex);

      friend istream & operator>> (istream &, complex&);
      friend ostream & operator<< (ostream &, complex);
};

// Inline funkcije

inline complex::complex(double r, double i) : re(r), im(i) {}

inline complex operator+ (complex c1, complex c2){
	return complex(c1.re+c2.re,c1.im+c2.im);
}

inline complex operator- (complex c1, complex c2){
	return complex(c1.re-c2.re,c1.im-c2.im);
}

inline complex operator* (complex c1, complex c2){
	return complex(c1.re*c2.re-c1.im*c2.im,c1.re*c2.im+c1.im*c2.re);
}

inline complex operator/ (complex c1, complex c2){
   double r=c2.re*c2.re+c2.im*c2.im;
	return complex((c1.re*c2.re+c1.im*c2.im)/r,(c2.re*c1.im-c2.im*c1.re)/r);
}

inline complex & complex::operator+= (complex c){
	re +=c.re; im +=c.im; return *this;
}

inline complex & complex::operator-= (complex c){
	re -=c.re; im -=c.im; return *this;
}

inline complex & complex::operator*= (complex c){
	return *this=*this*c;
}

inline complex & complex::operator/= (complex c){
	return *this=*this/c;
}

inline complex & complex::operator= (complex c) { re=c.re; im=c.im; return *this; }

  
inline complex operator+ (complex c) { return c; }

inline complex operator- (complex c) { return complex(-c.re,-c.im); }

inline int operator== (complex c1, complex c2){
	return (c1.re==c2.re) && (c1.im==c2.im);
}

inline int operator!= (complex c1, complex c2){
	return !(c1==c2);
}

inline double real (complex c) { return c.re; }

inline double imag (complex c) { return c.im; }

inline complex conj (complex c) { return complex(c.re,-c.im); }

inline double arg (complex c) { return c==complex(0,0) ? 0 : atan2(c.im,c.re); }

inline double norm (complex c) { return (c.re*c.re+c.im*c.im); }

inline double abs (complex c) { return sqrt(c.re*c.re+c.im*c.im); }

inline complex polar (double ro, double phi) { return complex(ro*cos(phi),ro*sin(phi)); }

inline complex exp (complex c) { return complex(exp(c.re)*cos(c.im),exp(c.re)*sin(c.im)); }

inline complex log (complex c) { return complex(log(abs(c)),arg(c)); }

inline complex log10 (complex c) { return log10(exp(1.0))*log(c); }

inline complex sqrt (complex c) {
   return complex(sqrt(abs(c))*cos(arg(c)/2.),sqrt(abs(c))*sin(arg(c)/2.));
}

// For power functions it was assumed : 0^0 == 1 ,  0^x=0 (when x !=0)

inline complex power (const complex a,int n){
   if(a==complex(0,0)) return n==0 ? complex(1,0) : complex(0,0);
   if(imag(a)==0) return real(a)<0 ? power(a,complex(n,0)) : complex(pow(real(a),n),0);
   return exp(complex((double)n,0)*log(a));
}

inline complex power (const complex a, double s) {
   if(a==complex(0,0)) return s==0 ? complex(1,0) : complex(0,0);
   if(imag(a)==0) return real(a)<0 ? power(a,complex(s,0)) : complex(pow(real(a),s),0);
	return exp(complex(s,0)*log(a));
}

inline complex power (double s, const complex a) {
   if(s==0) return a==complex(0,0) ? complex(1,0) : complex(0,0);
   if(s<0) return power(complex(s,0),a);
   if(imag(a)==0) return complex(pow(s,real(a)),0);
	return exp(a*complex(log(s),0));
}

inline complex power (const complex a1, const complex a2){
	double r1,th1,u2,v2,ro,phi;

   if(a1==complex(0,0)) return a2==complex(0,0) ? complex(1,0) : complex(0,0);

   r1=abs(a1); th1=arg(a1);
   u2=real(a2); v2=imag(a2);
   ro=pow(r1,u2)*exp(-v2*th1);
   phi=v2*log(r1)+u2*th1;
   return complex(ro*cos(phi),ro*sin(phi));
}

inline complex sin (complex a){
	return complex(sin(real(a))*cosh(imag(a)),cos(real(a))*sinh(imag(a)));
}

inline complex cos (complex a){
	return complex(cos(real(a))*cosh(imag(a)),-sin(real(a))*sinh(imag(a)));
}

inline complex tan (complex a) { return sin(a)/cos(a); }

inline complex sinh (complex a) {
	return complex(sinh(real(a))*cos(imag(a)),cosh(real(a))*sin(imag(a)));
}

inline complex cosh (complex a) {
	return complex(cosh(real(a))*cos(imag(a)),sinh(real(a))*sin(imag(a)));
}

inline complex tanh (complex a) { return sinh(a)/cosh(a); }

inline complex asin (complex a){
	const complex i(0,1);
   return i*log(i*a+sqrt(complex(1,0)-a*a));
}

inline complex acos (complex a){
	const complex i(0,1);
   return -i*log(a+i*sqrt(complex(1,0)-a*a));
}

inline complex atan (complex a){
	const complex i(0,1);
   return (i/complex(2,0))*log((i+a)/(i-a));;
}


inline istream & operator>> (istream &is, complex &c) {  return is >> c.re >> c.im; }

inline ostream & operator<< (ostream &os, complex c) {  return os << "(" << c.re << "," << c.im << ")"; }

#endif
