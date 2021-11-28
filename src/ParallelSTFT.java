import java.util.ArrayList;
import java.util.Arrays;
import java.util.concurrent.*;
import java.util.concurrent.Executors;

public class ParallelSTFT implements Callable {

    double[] x;
    double[] y;

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

    public static ArrayList<ArrayList<double[]>> parallelSTFT(double[] real, double[] imag, int periods, int numThreads) {

        int index = 0;
        int windowSize = real.length / periods;

        ArrayList<double[]> reals = new ArrayList<>();
        ArrayList<double[]> imags = new ArrayList<>();

        ArrayList<Future<ArrayList<double[]>>> res = new ArrayList<>();
        ExecutorService exec = Executors.newFixedThreadPool(numThreads);

        try {
            for(int i = 0; i < periods; i++) {
                Future<ArrayList<double[]>> f = exec.submit(new ParallelSTFT(Arrays.copyOfRange(real, index, index + windowSize),
                        Arrays.copyOfRange(imag, index, index + windowSize)));
                res.add(f);
                index += windowSize;
            }

            for(int i = 0; i < periods; i++) {
                ArrayList<double[]> r = res.get(i).get();

                reals.add(r.get(0));
                imags.add(r.get(1));
            }
        } catch(Exception e) {
            System.err.println(e);
        }

        exec.shutdown();

        ArrayList<ArrayList<double[]>> results = new ArrayList<>();
        results.add(reals);
        results.add(imags);
        return results;
    }

    public static void printResultsConsole(ArrayList<ArrayList<double[]>> res, int numWindows) {
        for(int i = 0; i < numWindows; i++) {
            System.out.println("Window: " + (i + 1));
            System.out.print("Reals: [ ");
            for(int j = 0; j < res.get(0).get(i).length; j++) {
                System.out.print(res.get(0).get(i)[j] + " ");
            }
            System.out.print("]\nImaginaries: [" );
            for(int j = 0; j < res.get(1).get(i).length; j++) {
                System.out.print(res.get(1).get(i)[j] + " ");
            }
            System.out.println("]");
        }
    }


    // Test functionality
    public static void main(String[] args) {
        int indx = 0;
        int numWindows = 8;

        int N = 1024;

        double[] re = new double[N];
        double[] im = new double[N];

        //Impulse
        for(int i = 0; i < N; i++) {
            if(i % 256 == 0) {
                re[i] = 2;
                im[i] = 0;
            } else if (i % 128 == 0) {
                re[i] = 1;
                im[i] = 0;
            } else
                re[i] = im[i] = 0;
        }

        ArrayList<ArrayList<double[]>> results = parallelSTFT(re, im, numWindows, 4);

        printResultsConsole(results, numWindows);

    }
}
