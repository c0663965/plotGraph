
// Fully written by Kihoon, Lee in 2017-08-01
// 향후 추가될 사항 : 선형회귀분석, 비중복난수, 중복원수 개수카운팅, txt화일 입출력, 그래프 그리기 

package javaapplication5;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

class MyMath {

    MyMath() {
    }

    public static int[] reverse(int[] x) {

        int temp;

        for (int i = 0; i < x.length / 2; i++) {
            temp = x[x.length - 1 - i];
            x[x.length - 1 - i] = x[i];
            x[i] = temp;
        }

        return x;
    }

    public static double[] reverse(double[] x) {

        double temp;

        for (int i = 0; i < x.length / 2; i++) {
            temp = x[x.length - 1 - i];
            x[x.length - 1 - i] = x[i];
            x[i] = temp;
        }

        return x;
    }

    public static Integer[] reverse(Integer[] x) {

        Integer temp;

        for (int i = 0; i < x.length / 2; i++) {
            temp = x[x.length - 1 - i];
            x[x.length - 1 - i] = x[i];
            x[i] = temp;
        }

        return x;
    }

    public static Double[] reverse(Double[] x) {

        Double temp;

        for (int i = 0; i < x.length / 2; i++) {
            temp = x[x.length - 1 - i];
            x[x.length - 1 - i] = x[i];
            x[i] = temp;
        }

        return x;
    }

    public static String[] reverse(String[] x) {

        String temp;

        for (int i = 0; i < x.length / 2; i++) {
            temp = x[x.length - 1 - i];
            x[x.length - 1 - i] = x[i];
            x[i] = temp;
        }

        return x;
    }

    public static int[][] plus(int[][]... args) {

        int n = args.length;

        int r = args[0].length;
        int c = args[0][0].length;

        int[][] A = newInt(r,c);
   
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < r; j++) {
                for (int k = 0; k < c; k++) {
                    A[j][k] += args[i][j][k];
                }
            }
        }

        return A;
    }
    
    public static double[][] plus(double[][]... args) {

        int n = args.length;

        int r = args[0].length;
        int c = args[0][0].length;

        double[][] A = newDouble(r,c);
   
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < r; j++) {
                for (int k = 0; k < c; k++) {
                    A[j][k] += args[i][j][k];
                }
            }
        }

        return A;
    }

    public static int[] plus(int[]... args) {

        int r = args.length;
        int c = args[0].length;

        int[] A = new int[args[0].length];
        Arrays.fill(A, 0);

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                A[j] += args[i][j];
            }
        }

        return A;
    }

    public static double[] plus(double[]... args) {

        int r = args.length;
        int c = args[0].length;

        double[] A = new double[args[0].length];
        Arrays.fill(A, 0);

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                A[j] += args[i][j];
            }
        }

        return A;
    }

    public static int[][] minus(int[][] A, int[][] B) {

        int r = A.length;
        int c = A[0].length;

        int[][] C = new int[r][c];

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                C[i][j] = A[i][j] - B[i][j];
            }
        }

        return C;
    }

    public static double[][] minus(double[][] A, double[][] B) {

        int r = A.length;
        int c = A[0].length;

        double[][] C = new double[r][c];

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                C[i][j] = A[i][j] - B[i][j];
            }
        }

        return C;
    }

    public static int[] minus(int[] A, int[] B) {

        int r = A.length;

        int[] C = new int[r];

        for (int i = 0; i < r; i++) {
            C[i] = A[i] - B[i];
        }

        return C;
    }

    public static double[] minus(double[] A, double[] B) {

        int r = A.length;

        double[] C = new double[r];

        for (int i = 0; i < r; i++) {
            C[i] = A[i] - B[i];
        }

        return C;
    }

     public static double[][] multiply(double[][]... args) { 
        int len = args.length;

        double[][][] temp = new double[len][][];
        temp[len-1] = new double[args[len-1].length][args[len-1][0].length];
        temp[len-1] = args[len-1];
      
        for (int i = len-1; i > 0; i--) {
            temp[i-1] = new double[args[i - 1].length][args[i][0].length];
            temp[i-1]=mul(args[i - 1], temp[i]);
        }

        return temp[0]; 
    }
     
    public static int[][] multiply(int[][]... args) {
        int len = args.length;

        int[][][] temp = new int[len][][];
        temp[len-1] = new int[args[len-1].length][args[len-1][0].length];
        temp[len-1] = args[len-1];
      
        for (int i = len-1; i > 0; i--) {
            temp[i-1] = new int[args[i - 1].length][args[i][0].length];
            temp[i-1]=mul(args[i - 1], temp[i]);
        }

        return temp[0]; 
    }
        
    private static double[][] mul(double[][] A, double[][] B) {

        int r = A.length;
        int c = B[0].length;
        int u = B.length;
        double[][] C = newDouble(r,c);

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                for (int k = 0; k < u; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return C;
    }
    
    private static int[][] mul(int[][] A, int[][] B) {

        int r = A.length;
        int c = B[0].length;
        int u = B.length;
        int[][] C = newInt(r,c);

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                for (int k = 0; k < u; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return C;
    }
 
    public static int[] multiply(int[][] A, int[] x) {

        int r = A.length;
        int c = A[0].length;

        int[] y = new int[r];
        Arrays.fill(y,0);
        
        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                y[i] += A[i][j] * x[j];
            }
        }

        return y;
    }

    public static double[] multiply(double[][] A, double[] x) {

        int r = A.length;
        int c = A[0].length;

        double[] y = new double[r];
        Arrays.fill(y,0);

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                y[i] += A[i][j] * x[j];
            }
        }

        return y;
    }

    public static int[][] multiply(int[] x, int[] y) {

        int r = x.length;
        int c = y.length;

        int[][] Z = new int[r][c];

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                Z[i][j] = x[i] * y[j];
            }
        }

        return Z;
    }

    public static double[][] multiply(double[] x, double[] y) {

        int r = x.length;
        int c = y.length;

        double[][] Z = new double[r][c];

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                Z[i][j] = x[i] * y[j];
            }
        }

        return Z;
    }

    public static int[] multiply(int k, int[] x) {

        int r = x.length;

        int[] kx = new int[r];

        for (int i = 0; i < r; i++) {
            kx[i] = k * x[i];
        }

        return kx;
    }

    public static double[] multiply(double k, double[] x) {

        int r = x.length;

        double[] kx = new double[r];

        for (int i = 0; i < r; i++) {
            kx[i] = k * x[i];
        }

        return kx;
    }

    public static Integer[] multiply(Integer k, Integer[] x) {

        int r = x.length;

        Integer[] kx = new Integer[r];

        for (int i = 0; i < r; i++) {
            kx[i] = k * x[i];
        }

        return kx;
    }

    public static Double[] multiply(Double k, Double[] x) {

        int r = x.length;

        Double[] kx = new Double[r];

        for (int i = 0; i < r; i++) {
            kx[i] = k * x[i];
        }

        return kx;
    }

    public static int[][] multiply(int k, int[][] A) {

        for (int[] r : A) {
            for (int j = 0; j < A[0].length; j++) {
                r[j] = k * r[j];
            }
        }
        return A;
    }

    public static double[][] multiply(double k, double[][] A) {

        for (double[] r : A) {
            for (int j = 0; j < A[0].length; j++) {
                r[j] = k * r[j];
            }
        }
        return A;
    }

    public static Integer[][] multiply(Integer k, Integer[][] A) {

        for (Integer[] r : A) {
            for (int j = 0; j < A[0].length; j++) {
                r[j] = k * r[j];
            }
        }
        return A;
    }

    public static Double[][] multiply(Double k, Double[][] A) {

        for (Double[] r : A) {
            for (int j = 0; j < A[0].length; j++) {
                r[j] = k * r[j];
            }
        }
        return A;
    }

    public static int product(int[] x, int[] y) {

        int sum = 0;

        for (int i = 0; i < x.length; i++) {
            sum += x[i] * y[i];
        }

        return sum;
    }

    public static double product(double[] x, double[] y) {

        double sum = 0;

        for (int i = 0; i < x.length; i++) {
            sum += x[i] * y[i];
        }

        return sum;
    }

    public static double[][] transpose(double[][] A) {
        int r = A.length;
        int c = A[0].length;

        double[][] T = new double[c][r];

        for (int i = 0; i < c; i++) {
            for (int j = 0; j < r; j++) {
                T[i][j] = A[j][i];
            }
        }

        return T;
    }

    public static int[][] transpose(int[][] A) {
        int r = A.length;
        int c = A[0].length;

        int[][] T = new int[c][r];

        for (int i = 0; i < c; i++) {
            for (int j = 0; j < r; j++) {
                T[i][j] = A[j][i];
            }
        }

        return T;
    }

    public static int sum(int[] arr) {
        int sum = 0;

        for (int e : arr) {
            sum += e;
        }

        return sum;
    }

    public static double avg(int[] arr) {
        return (double) sum(arr) / arr.length;
    }
    
    public static double stdev(int [] arr){
        
        int n=arr.length;
        double average=avg(arr);
            
        double temp=0;
        for (int i=0;i<n;i++)
            temp+=Math.pow((arr[i]-average),2);
        
        return Math.sqrt(temp/n); 
    }
 
    public static double sum(double[] arr) {
        double sum = 0;

        for (double e : arr) {
            sum += e;
        }

        return sum;
    }

    public static double avg(double[] arr) {
        return (double) sum(arr) / arr.length;
    }
    
    public static double stdev(double[] arr){
        
        int len=arr.length;
        double average=avg(arr);
            
        double temp=0;
        for (int i=0;i<len;i++)
            temp+=Math.pow((arr[i]-average),2);
        
        return Math.sqrt(temp/len); 
    }
    
    public static Integer sum(Integer[] arr) {
        Integer sum = 0;

        for (Integer e : arr) {
            sum += e;
        }

        return sum;
    }

    public static Integer avg(Integer[] arr) {
        return sum(arr) / arr.length;
    }
    
    public static double stdev(Integer[] arr){
        
        int len=arr.length;
        Integer average=avg(arr);
            
        double temp=0;
        for (int i=0;i<len;i++)
            temp+=Math.pow((arr[i]-average),2);
        
        return Math.sqrt(temp/len); 
    }

    public static Double sum(Double[] arr) {
        Double sum = 0.0;

        for (Double e : arr) {
            sum += e;
        }

        return sum;
    }

    public static Double avg(Double[] arr) {
        return sum(arr) / arr.length;
    }
    
    public static Double stdev(Double[] arr){
        
        int len=arr.length;
        Double average=avg(arr);
            
        Double temp=0.0;
        for (int i=0;i<len;i++)
            temp+=Math.pow((arr[i]-average),2);
        
        return Math.sqrt(temp/len); 
    }

    public static int max(int[] arr) {

        Arrays.sort(arr);

        return arr[arr.length - 1];
    }

    public static Integer max(Integer[] arr) {

        Arrays.sort(arr);

        return arr[arr.length - 1];
    }

    public static double max(double[] arr) {

        Arrays.sort(arr);

        return arr[arr.length - 1];
    }

    public static Double max(Double[] arr) {

        Arrays.sort(arr);

        return arr[arr.length - 1];
    }

    public static int min(int[] arr) {

        Arrays.sort(arr);

        return arr[0];
    }

    public static Integer min(Integer[] arr) {

        Arrays.sort(arr);

        return arr[0];
    }

    public static double min(double[] arr) {

        Arrays.sort(arr);

        return arr[0];
    }

    public static Double min(Double[] arr) {

        Arrays.sort(arr);

        return arr[0];
    }

    public static int[][] newInt(int r, int c) {
        int[][] A = new int[r][c];

        for (int p[] : A) {
            Arrays.fill(p, 0);
        }

        return A;
    }

    public static double[][] newDouble(int r, int c) {
        double[][] A = new double[r][c];

        for (double p[] : A) {
            Arrays.fill(p, 0);
        }

        return A;
    }

    public static int[] newInt(int n) {
        int[] A = new int[n];
        Arrays.fill(A, 0);

        return A;
    }

    public static double[] newDouble(int n) {
        double[] A = new double[n];
        Arrays.fill(A, 0.0);

        return A;
    }

    public static Integer[] newInteger(int n) {
        Integer[] A = new Integer[n];
        Arrays.fill(A, 0);

        return A;
    }

    public static Double[] newDOUBLE(int n) {
        Double[] A = new Double[n];
        Arrays.fill(A, 0.0);

        return A;
    }

    public static void show(int[][] a) {
        for (int[] r : a) {
            System.out.println("");
            for (int j = 0; j < r.length; j++) {
                System.out.print(r[j] + "  \t");
            }
        }

        System.out.println("");
    }

    public static void show(Integer[][] a) {
        for (Integer[] r : a) {
            System.out.println("");
            for (Integer r1 : r) {
                System.out.print(r1 + " ");
            }
        }

        System.out.println("");
    }

    public static void show(double[][] a) {

        NumberFormat formatter = new DecimalFormat("#0.00");

        for (double[] r : a) {
            System.out.println("");
            for (int j = 0; j < r.length; j++) {
                System.out.print(formatter.format(r[j]) + " ");
            }
        }

        System.out.println("");
    }

    public static void show(Double[][] a) {

        NumberFormat formatter = new DecimalFormat("#0.00");

        for (Double[] r : a) {
            System.out.println("");
            for (Double r1 : r) {
                System.out.print(formatter.format(r1) + " ");
            }
        }

        System.out.println("");
    }

    public static void show(int[] m) {
        System.out.println(Arrays.toString(m));
    }

    public static void show(Integer[] m) {
       System.out.println(Arrays.toString(m));
    }

    public static void show(double[] m) {
       System.out.println(Arrays.toString(m));
    }

    public static void show(Double[] m) {
        System.out.println(Arrays.toString(m));
    }

    public static void show(List<List> a) {
        for (int i = 0; i < a.size(); i++) {
            System.out.println();
            for (int j = 0; j < a.get(i).size(); j++) {
                System.out.print(a.get(i).get(j) + "   \t");
            }
        }

        System.out.println("");
    }

    public static int[] randInt(int len) {
        int[] arr = newInt(len);

        Random rd = new Random();

        for (int i = 0; i < len; i++) {
            arr[i] = rd.nextInt();
        }

        return arr;
    }
    
    public static double[] randDouble(int len) {
        double[] arr = newDouble(len);

        Random rd = new Random();

        for (int i = 0; i < len; i++) {
            arr[i] = rd.nextDouble();
        }

        return arr;
    }
    
    public static int[] randInt(int len, int min, int max) {
        int[] rn = newInt(len);

        Random rd = new Random();

        for (int i = 0; i < len; i++) {
            rn[i] = rd.nextInt(max - min + 1) + min;
        }

        return rn;
    }

    public static double[] randDouble(int len, int min, int max) {
        double[] rn = newDouble(len);

        Random rd = new Random();

        for (int i = 0; i < len; i++) {
            rn[i] = rd.nextDouble() * (max - min) + min;
        }

        return rn;
    }

    public static double[][] randDouble(int r, int c) {

        double[][] mat = newDouble(r, c);

        Random rd = new Random();

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                mat [i][j] = rd.nextDouble();
            }
        }

        return mat ;
    }

    public static double[][] randDouble(int r, int c, int min, int max) {

        double[][] mat = newDouble(r, c);

        Random rd = new Random();

        for (int i = 0; i < r; i++) {
            for (int j = 0; j < c; j++) {
                mat[i][j] = rd.nextDouble() * (max - min) + min;
            }
        }

        return mat;
    }

    public static int[] naturalNumber(int max) {  //max까지의 자연수

        int[] arr = newInt(max);

        for (int i = 0; i < max; i++) {
            arr[i] = i + 1;
        }

        return arr;
    }
    
    public static int[] shuffle(int[] arr) { 

        Random rd = new Random();

        int max=arr.length;
        
        int n = 0;
        int temp, index1, index2;

        while (n++ < max) {
            index1 = rd.nextInt(arr.length);
            index2 = rd.nextInt(arr.length);

            temp = arr[index1];
            arr[index1] = arr[index2];
            arr[index2] = temp;
        }

        return arr;
    }
    
    public static double[] shuffle(double[] arr) { 

        Random rd = new Random();

        int max=arr.length;
        
        int n = 0;
        int index1, index2;
        double temp;
        
        while (n++ < max) {
            index1 = rd.nextInt(arr.length);
            index2 = rd.nextInt(arr.length);

            temp = arr[index1];
            arr[index1] = arr[index2];
            arr[index2] = temp;
        }

        return arr;
    }

    public static Integer[] int2Integer(int[] arr) {

        Integer[] r = new Integer[arr.length];

        for (int i = 0; i < arr.length; i++) {
            r[i] = arr[i];
        }

        return r;
    }

    public static int[] Integer2int(Integer[] arr) {

        int[] r = new int[arr.length];

        for (int i = 0; i < arr.length; i++) {
            r[i] = arr[i];
        }

        return r;
    }

    public static Double[] double2Double(double[] arr) {

        Double[] r = new Double[arr.length];

        for (int i = 0; i < arr.length; i++) {
            r[i] = arr[i];
        }

        return r;
    }

    public static double[] Double2double(Double[] arr) {

        double[] r = new double[arr.length];

        for (int i = 0; i < arr.length; i++) {
            r[i] = arr[i];
        }

        return r;
    }

    public static double[] linespan(double start, double finale, double step) {
        List<Double> temp = new ArrayList();

        int i = 0;

        temp.add(start);

        while (temp.get(i++) < finale) {
            temp.add(start + i * step);
        }

        double[] t = newDouble(temp.size());

        int j = 0;

        for (Double r : temp) {
            t[j++] = r;
        }

        return t;
    }

    public static double[][] inv(double[][] A) {

        int row = A.length;
        int col = 2 * row;
        
        double[][] T = new double[row][col];
        double temp;

        for (int i = 0; i < row; i++) {
            System.arraycopy(A[i], 0, T[i], 0, row);

            for (int j = row; j < col; j++) {
                if (j == row + i) {
                    T[i][j] = 1;
                } else {
                    T[i][j] = 0;
                }
            }
        }
        
        /**
         * ******************* Downward*************************
         */
        for (int k = 0; k < row - 1; k++) {
            temp = T[k][k];
            if (temp == 0) {
                T = mix(T, k);
                temp = T[k][k];
            }
            for (int j = 0; j < col; j++) {
                T[k][j] = T[k][j] / temp;
            }
            for (int i = k + 1; i < row; i++) {
                temp = T[i][k];
                for (int j = 0; j < col; j++) {
                    T[i][j] = -T[i][j] + T[k][j] * temp;
                }
            }
        }
        /**
         * ******************* Upward****************************
         */
        for (int k = row - 1; k >= 0; k--) {
            temp = T[k][k];

            if (temp == 0) {
                T = mix(T, k);
                temp = T[k][k];
            }

            for (int j = 0; j < col; j++) {
                T[k][j] = T[k][j] / temp;
            }
            for (int i = k - 1; i >= 0; i--) {
                temp = T[i][k];
                for (int j = 0; j < col; j++) {
                    T[i][j] = -T[i][j] + T[k][j] * temp;
                }
            }
        }

        double[][] inv_A = new double[row][row];

        for (int i = 0; i < row; i++) {
            for (int j = 0; j < row; j++) {
                inv_A[i][j] = T[i][j + row];
            }
        }

        return inv_A;
    }

    private static double[][] mix(double a[][], int k){ //To avoid the case of a[k][k]==0
   
        int row = a.length;
        int col = a[0].length;
        
        double temp = 0;

        for (int j = 0; j < col; j++) {
            for (int i = 0; i < row; i++) {
                temp += a[i][j];
            }
            a[k][j] = temp;
            temp = 0;
        }
        
        return a;
    }
    
    public static double det(double a[][]) {

        double det;
        det = 0;
        if (a.length == 2) {
            det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        } else {
            for (int j = 0; j < a.length; j++) {
                det = det + a[0][j] * (double) Math.pow(-1, 0 + j) * det(cofactor(a, 0, j));
            }
        }   // Recursive calling of the functino Det.
        return det;
    }
    
    private static double[][] cofactor(double a[][], int i, int j) {

        double temp1;

        double[][] b = new double[a.length][a.length];

        int n=0;
        
        for (double[] r:a)
            b[n++] = Arrays.copyOf(r,a.length);

        while (j > 0) {
            for (int k = 0; k < a.length; k++) {
                temp1 = b[k][j];
                b[k][j] = b[k][j - 1];
                b[k][j - 1] = temp1;
            }

            j--;
        }
        while (i > 0) {
            for (int k = 0; k < a.length; k++) {
                temp1 = b[i][k];
                b[i][k] = b[i - 1][k];
                b[i - 1][k] = temp1;
            }
            i--;
        }

        /* Construction of the new smaller confactor matrix */
        double[][] cof = new double[a.length - 1][a.length - 1];

        for (int p = 0; p < a.length - 1; p++) {
            for (int q = 0; q < a.length - 1; q++) {
                cof[p][q] = b[p + 1][q + 1];
            }
        }

        return cof;
    }
    
    public static int det(int a[][]) {

        int det;
        det = 0;
        if (a.length == 2) {
            det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        } else {
            for (int j = 0; j < a.length; j++) {
                det = det + a[0][j] * (int)Math.pow(-1,j) * det(cofactor(a, 0, j));
            }
        }   // Recursive calling of the functino Det.
        return det;
    }
    
    private static int[][] cofactor(int a[][], int i, int j) {

        int temp1;

        int[][] b = new int[a.length][a.length];

        int n=0;
        
        for (int[] r:a)
            b[n++] = Arrays.copyOf(r,a.length);

        while (j > 0) {
            for (int k = 0; k < a.length; k++) {
                temp1 = b[k][j];
                b[k][j] = b[k][j - 1];
                b[k][j - 1] = temp1;
            }

            j--;
        }
        
        while (i > 0) {
            for (int k = 0; k < a.length; k++) {
                temp1 = b[i][k];
                b[i][k] = b[i - 1][k];
                b[i - 1][k] = temp1;
            }
            i--;
        }

        /* Construction of the new smaller confactor matrix */
        int[][] cof = new int[a.length - 1][a.length - 1];

        for (int p = 0; p < a.length - 1; p++) {
            for (int q = 0; q < a.length - 1; q++) {
                cof[p][q] = b[p + 1][q + 1];
            }
        }

        return cof;
    }
    
    public static List<List<Double>> RungeKutta4(double[][] A, double[] x0, double t0, double tf, double dt ){
        List<List<Double>> data = new ArrayList();
        
        double[] k = newDouble(x0.length);
        double[] h = {0, dt/2.0, dt/2.0, dt};
        double[] c = {dt/6.0, dt/3.0, dt/3.0, dt/6.0};
        
        double[] x = newDouble(x0.length);
        double[] dx = newDouble(x0.length);
        double[] temp = newDouble(1+x0.length); 
        
        System.arraycopy(x0, 0, temp, 1, x0.length);
       
        data.add(Arrays.asList(double2Double(temp)));
        
        while(t0<tf)
        {    
            for (int i = 0 ; i<4 ; i++)
            {
               x=plus(x0,multiply(h[i],k)); //dt*time 
               
               k=multiply(A,x);   //dt
               
               dx=plus(dx,multiply(c[i],k)); //dt*time         
            }
            
            x0=plus(x0,dx);
            
            Arrays.fill(dx,0);
            
            t0+=dt;
            
            temp[0]=t0;
            
            System.arraycopy(x0, 0, temp, 1, x0.length);
          
            data.add(Arrays.asList(double2Double(temp)));
        }
        
        return data;
    }
    
    public static double[] concat(double[]...args) {
        
        int len = args.length;
        int sum=0;
        
        for (double r[]:args)
            sum+=r.length;
    
        double[] v = newDouble(sum);
        
        for (int i=0,desPos=0;i<len;i++){
            System.arraycopy(args[i],0,v,desPos,args[i].length);
            desPos+=args[i].length;
        }
                    
        return v;
    }
    
    public static int[] concat(int[]...args) {
        
        int len = args.length;
        int sum=0;
  
        for (int r[]:args)
            sum+=r.length;
        
        int[] v = newInt(sum);
        
        for (int i=0,desPos=0;i<len;i++){
            System.arraycopy(args[i],0,v,desPos,args[i].length);
            desPos+=args[i].length;
        }
                    
        return v;
    }
    
    public static double[] concat(double n, double[] v) {
        double[] w = newDouble(1+v.length);
        w[0]=n;
        System.arraycopy(v, 0, w, 1, v.length);
        return w;
    }
    
    public static double[] concat(double[] v, double n) {
        double[] w = newDouble(v.length+1);
        System.arraycopy(v, 0, w, 0, v.length);
        w[w.length-1]=n;
   
        return w;
    }
    
    public static int[] concat(int n, int[] v) {
        int[] c = newInt(1+v.length);
        c[0]=n;
        System.arraycopy(v, 0, c, 1, v.length);
        
        return c;
    }
    
    public static int[] concat(int[] v, int n) {
        int[] w = new int[v.length+1];
        System.arraycopy(v, 0, w, 0, v.length);
        w[w.length-1]=n;

        return w;
    }
    
    public static int[][] reshape(int[] v, int...args){ //reshape(new int[]{1,2,3,4,5,6,7,8,9},3); 이런식으로도 사용가능함
        
        int len=args.length;
        int r,c;
        
        if (len==1){
            r=args[0];
            c=v.length/r;
        }else{
            r=args[0];
            c=args[1];
        }
        
        int[][] a = new int[r][c];
        
        for (int i=0;i<r;i++)
            for (int j=0;j<c;j++)
                a[i][j]=v[c*i+j];
        
        return a;
    }
    
    public static double[][] reshape(double[] v, int...args){
        
        int len=args.length;
        int r,c;
        
        if (len==1){
            r=args[0];
            c=v.length/r;
        }else{
            r=args[0];
            c=args[1];
        }
        
        double[][] a = new double[r][c];
        
        for (int i=0;i<r;i++)
            for (int j=0;j<c;j++)
                a[i][j]=v[c*i+j];
        
        return a;
    }
    
    
    
    public static int[][] reshape(int[][] arr, int...args){
       
        int len=args.length;
        int r,c;
        
        if (len==1){
            r=args[0];
            c=arr.length*arr[0].length/r;
        }else{
            r=args[0];
            c=args[1];
        }
    
        int size = arr.length*arr[0].length;
        int[] v = new int[size];

        int i=0;

        for (int[] row : arr)
           for (int col :row)
                v[i++]=col;

        int[][] a = new int[r][c];

        for (int j=0;j<r;j++)
            for (int k=0;k<c;k++)
                a[j][k]=v[c*j+k];

        return a;
    }
    
    public static double[][] reshape(double[][] arr, int...args){
        
        int len=args.length;
        int r,c;
        
        if (len==1){
            r=args[0];
            c=arr.length*arr[0].length/r;
        }else{
            r=args[0];
            c=args[1];
        }
    
        int size = arr.length*arr[0].length;
        double[] v = new double[size];

        int i=0;

        for (double[] row : arr)
           for (double col :row)
                v[i++]=col;

        double[][] a = new double[r][c];

        for (int j=0;j<r;j++)
            for (int k=0;k<c;k++)
                a[j][k]=v[c*j+k];

        return a;
    }
    
    public static void plot(List<Double> x, List<Double> y){ //jfreechart-1.0.19.jar and jcommon-1.023.jar
            
        final XYSeries data = new XYSeries("");

        for (int i=0;i<x.size();i++)
            data.add(x.get(i),y.get(i));

        final XYSeriesCollection dataset = new XYSeriesCollection( );

        dataset.addSeries(data);

        JFreeChart plot = ChartFactory.createXYLineChart(
            "", 
            "",
            "", 
            dataset,
            PlotOrientation.VERTICAL, 
            false, false, false);

        ChartFrame frame = new ChartFrame("Results",plot);
        frame.pack();
        frame.setVisible(true);
    }
    
    public static void plot(List<List<Double>> data){
        final XYSeriesCollection dataset = new XYSeriesCollection( );
        
        for(int i=1;i<data.get(0).size();i++){
            final XYSeries temp= new XYSeries("Legend"+i);
            for (List<Double> r : data)
                temp.add(r.get(0),r.get(i));
            dataset.addSeries(temp);
        }
        
        JFreeChart plot = ChartFactory.createXYLineChart(
            "Linear System Response", 
            "time(s)",
            "State", 
            dataset,
            PlotOrientation.VERTICAL, 
            true, true, false);

        ChartFrame frame = new ChartFrame("Results", plot);
        frame.pack();
        frame.setVisible(true);
    }
    
    public static void plot(double[] x, double[] y){
        List<Double> X= Arrays.asList(double2Double(x));
        List<Double> Y= Arrays.asList(double2Double(y));
        plot(X,Y);
    }
    
    public static void plot(double[][] data){
        List<List<Double>> temp = new ArrayList();
        for (double r[]:data)
            temp.add(Arrays.asList(double2Double(r)));
       plot(temp);
    }
}