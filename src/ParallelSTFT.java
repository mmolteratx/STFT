import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.*;
import java.util.concurrent.Executors;

public class ParallelSTFT implements Callable {

    double[] x;
    double[] y;
    ArrayList<ArrayList<double[]>> results;

    public ParallelSTFT(double[] x, double[] y) {
        this.x = x;
        this.y = y;
    }

    public ArrayList<double[]> call() {
        FFT fft = new FFT(x.length);
        fft.fft(x, y);

        ArrayList<double[]> res = new ArrayList<>();
        res.add(x);
        res.add(y);

        return res;
    }

    // Test functionality
    // TODO: Combine results from all threads
    public static void main(String[] args) {
        int indx = 0;
        int numWindows = 3;

        int N = 96;

        double[] re = new double[N];
        double[] im = new double[N];

        //Impulse
        for(int i = 0; i < N; i++) {
            if(i % 32 == 0) {
                re[i] = 1;
                im[i] = 0;
            } else
                re[i] = im[i] = 0;
        }
        int indx2 = 0;

        ArrayList<Future<ArrayList<double[]>>> res = new ArrayList<>();
        ExecutorService exec = Executors.newFixedThreadPool(4);

        try {
            while(indx < 3) {
                Future<ArrayList<double[]>> f = exec.submit(new ParallelSTFT(Arrays.copyOfRange(re, indx2, indx2 + 32),
                        Arrays.copyOfRange(im, indx2, indx2 + 32)));
                res.add(f);
                indx2 += 32;
                indx++;
            }

            indx = 0;

            while(indx < 3) {
                ArrayList<double[]> r = res.get(indx).get();
                System.out.println(r.get(0)[0]);
                System.out.println(r.get(1)[0]);
                indx++;
            }
        } catch(Exception e) {
            System.err.println(e);
        }




    }
}
