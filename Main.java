package com.lr1;
import java.lang.Math;

import static java.lang.Math.pow;
import static java.lang.Math.sqrt;

public class Main {

    public static double M1 = 4.60517018598809;
    public static double M2 = 19.2075924419136;
    public static double M4 = 195.270867879347;
    //public static double M8 = 331862.811526352983;
    public static double M10 = 1.39157296815225*pow(10,7);
    public static double m10 = 567262.155839182625;
    public static double EPS = pow(10,(-5));
    public static double EXACT_INTEGRAL = -0.1978168271761324;


    // МЕТОД ПРАВЫХ ПРЯМОУГОЛЬНИКОВ //
    public static double get_xk_rightR(double h, double k, double a){
        double xk = a + k*h;
        return ((pow(xk,2)-1)*pow(10,(-2*xk)));
    }
    public static double get_N_rightR(double a, double b){
        double N = pow(b-a,2)*M1/(2*EPS);
        return (int)(N+1);
    }
    public static double get_h_rightR(double N, double a, double b){
        return (b-a)/N;
    }
    public static double method_rightR(double h, double N, double a) {
        double sum = 0;
        for (int i=0; i<=N-1; i++){
            sum = sum + get_xk_rightR(h,i,a);
        }
        return h*sum;
    }
    public static double residuals_rightR(double N, double a, double b){
        return pow(b-a,2)*M1/(2*N);
    }
    public static double Runge_rule(double I1, double I2, double h1, double h2){
        int m = 1;
        return (I2-I1)/(1-pow(h2/h1,m));
    }
    public static double exact_integral_Runge(double I1, double I2, double h1, double h2){
        double res = Runge_rule(I1, I2, h1, h2);
        return I1+res;
    }


    // МЕТОД СРЕДНИХ ПРЯМОУГОЛЬНИКОВ //
    public static double get_xk_middleR(double h, double k, double a){
        double xk = a + k*h + h/2;
        return ((pow(xk,2)-1)*pow(10,(-2*xk)));
    }
    public static double get_N_middleR(double a, double b){
        double N = sqrt(pow(b-a,3)*M2/(24*EPS));
        return (int)(N+1);
    }
    public static double get_h_middleR(double N, double a, double b){
        return (b-a)/N;
    }
    public static double method_middleR(double h, double a, double b) {
        double N = get_N_middleR(a,b);
        double sum = 0;
        for (int i=0; i<=N-1; i++){
            sum = sum + get_xk_middleR(h,i,a);
        }
        return h*sum;
    }
    public static double residuals_middleR(double N, double a, double b){
        return pow(b-a,3)*M2/(24*pow(N,2));
    }


    // МЕТОД СИМПСОНА //
    public static double get_xk_Simpson(double h, double k, double a){
        double xk = a + k*h;
        return ((pow(xk,2)-1)*pow(10,(-2*xk)));
    }
    public static double get_N_Simpson(double a, double b){
        double N = pow((pow(b-a,5)*M4/(2880*EPS)), 0.25);
        return 2*((int)N+1);
    }
    public static double get_h_Simpson(double N, double a, double b){
        return (b-a)/N;
    }
    public static double method_Simpson(double h, double a, double b){
        double sum = 0;
        double N = get_N_Simpson(a,b);
        double x0 = get_xk_Simpson(h,0, a);
        double xN = get_xk_Simpson(h, N, a);
        double sumOdd = 0;
        double sumEven = 0;
        for(int i=1; i<=N-1; i++){
            if(i%2 == 0){
                sumEven = sumEven + get_xk_Simpson(h,i,a);
            }
            else{
                sumOdd = sumOdd + get_xk_Simpson(h,i,a);
            }
        }
        sum = x0 + xN + 4*sumOdd +2*sumEven;
        return h*sum/3;
    }
    public static double residuals_Simpson(double N, double a, double b){
        return M4/(180*pow(N,4));
    }

    // КФ НАСТ ГАУССА //
    public static int factorial (int n){
        if (n == 0)
            return 1;
        else
            return (n * factorial(n-1));
    }
    public static double new_x (double t, double a, double b){
        double x = (b-a)*t/2 + (b+a)/2;
        return x;
    }
    public static double get_tk_Gauss(double t){
        return ((pow(t,2)-1)*pow(10,(-2*t)));
    }
    public static double KF_NAST_Gauss3(double a, double b, double n){  // для 3 узлов
        double t1 = 0.86113631;
        double t2 = 0.33998104;
        double t3 = -0.33998104;
        double t4 = -0.86113631;
        double A1 = 0.34785484;
        double A2 = 0.65214516;
        double A3 = 0.65214516;
        double A4 = 0.34785484;
        double integral = ((b-a)/2)*(A1*get_tk_Gauss(new_x(t1,a,b)) + A2*get_tk_Gauss(new_x(t2,a,b)) +
                                     A3*get_tk_Gauss(new_x(t3,a,b)) + A4*get_tk_Gauss(new_x(t4,a,b)));
        return integral;
    }
    public static double KF_NAST_Gauss4(double a, double b, double n){  // чит для 4 узлов
        double t0 = 0;
        double t1 = 0.90617985;
        double t2 = 0.53846931;
        double t3 = -0.53846931;
        double t4 = -0.90617985;
        double A0 = 0.56888889;
        double A1 = 0.23692688;
        double A2 = 0.47862868;
        double A3 = 0.47862868;
        double A4 = 0.23692688;
        double integral = ((b-a)/2)*(A1*get_tk_Gauss(new_x(t1,a,b)) + A2*get_tk_Gauss(new_x(t2,a,b)) +
                A3*get_tk_Gauss(new_x(t3,a,b)) + A4*get_tk_Gauss(new_x(t4,a,b)) + A0*get_tk_Gauss(new_x(t0,a,b)));
        return integral;
    }
    public static double residuals_Gauss(int n,double a, double b){
        return M10*pow(b-a,2*n+3)*pow(factorial(n+1), 4)/((2*n+3)*pow(factorial(2*n+2), 3));
    }

    public static void main(String[] args) {

        double a = 0;
        double b = 1;

        // МЕТОД ПРАВЫХ ПРЯМОУГОЛЬНИКОВ //
        System.out.println("Правило Рунге");
        double N1_rR = get_N_rightR(a,b);
        double N2_rR = N1_rR*2;
        double h1_rR = get_h_rightR(N1_rR,a,b);
        double h2_rR = h1_rR/2;
        double integral1_rightR = method_rightR(h1_rR, N1_rR, a);
        double integral2_rightR = method_rightR(h2_rR, N2_rR, a);
        System.out.println("N1 = " + N1_rR);
        System.out.println("N2 = " + N2_rR);
        System.out.println("h1 = " + h1_rR);
        System.out.println("h2 = " + h2_rR);
        System.out.println("Приближенное значение интеграла по h1: " + integral1_rightR);
        System.out.println("Приближенное значение интеграла по h2: " + integral2_rightR);
        double res_Runge = Runge_rule(integral1_rightR, integral2_rightR, h1_rR, h2_rR);
        System.out.println("Точность: " + res_Runge);
        //double exact1 = exact_integral_Runge(integral1_rightR, integral2_rightR, h1_rR, h2_rR);
        //double exact2 = exact_integral_Runge(integral2_rightR, integral1_rightR, h1_rR, h2_rR);
        //System.out.println("Точное значение интеграла по h1: " + exact1);
        //System.out.println("Точное значение интеграла по h2: " + exact2);
        System.out.println();

        // МЕТОД СРЕДНИХ ПРЯМОУГОЛЬНИКОВ //
        System.out.println("Метод средних прямоугольников");
        double N_mR = get_N_middleR(a,b);
        double h_mR = get_h_middleR(N_mR,a,b);
        double integral_middleR = method_middleR(h_mR,a,b);
        System.out.println("N = " + N_mR);
        System.out.println("h = " + h_mR);
        System.out.println("Приближенное значение интеграла: " + integral_middleR);
        double res_middleR = residuals_middleR(N_mR,a,b);
        System.out.println("Точность: " + res_middleR);
        System.out.println();

        // МЕТОД СИМПСОНА //
        System.out.println("Метод Симпсона");
        double N_S = get_N_Simpson(a,b);
        double h_S = get_h_Simpson(N_S,a,b);
        double integral_Simpson = method_Simpson(h_S,a,b);
        System.out.println("N = " + N_S);
        System.out.println("h = " + h_S);
        System.out.println("Приближенное значение интеграла: " + integral_Simpson);
        double res_Simpson = residuals_Simpson(N_S,a,b);
        System.out.println("Точность: " + res_Simpson);
        System.out.println();

        // КФ НАСТ Гаусса //
        System.out.println("Квадратурная формула НАСТ типа Гаусса");
        double integral_Gauss = KF_NAST_Gauss4(a,b, 4);
        System.out.println("Приближенное значение интеграла: " + integral_Gauss);
        double res_Gauss = residuals_Gauss(4, a, b);
        System.out.println("Точность: " + res_Gauss);
    }
}