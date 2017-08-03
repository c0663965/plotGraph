package javaapplication5;

import java.util.ArrayList;
import java.util.List;
import static javaapplication5.MyMath.*;


public class XYLineChart_image {

   public static void main( String[ ] args )throws Exception {
       
        double[][] A={{-0.3,-0.4},{1,0}};
        double x[] = {1,1};

        double t0=0;
        double tf=50;
        double dt=0.1;

        List<List<Double>> data = RungeKutta4(A, x, t0, tf, dt);
       
        List<Double> t = new ArrayList();
        List<Double> v = new ArrayList();
        List<Double> s = new ArrayList();
        
        for (List<Double> r : data){
            t.add(r.get(0));
               v.add(r.get(1));
                    s.add(r.get(2));
        }
       plot(data);
   }
}