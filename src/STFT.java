import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class STFT {
    int numWindows;

    public STFT(int numWindows) {
        this.numWindows = numWindows;
    }

    public ArrayList<ArrayList<double[]>> stft(double[] x, double[] y) {
        int size = x.length;
        int windowSize = size / numWindows;
        int m = (int)(Math.log(windowSize) / Math.log(2));
        ArrayList<double[]> realRes = new ArrayList<>();
        ArrayList<double[]> imagRes = new ArrayList<>();

        if(windowSize != (1<<m)) {
            throw new RuntimeException("FFT length must be power of 2");
        }

        int indx = 0;
        int currentWindow = 0;

        FFT fft = new FFT(windowSize);

        while(currentWindow < numWindows) {
            System.out.println(indx + " " + windowSize);
            double[] real = Arrays.copyOfRange(x, indx, indx + windowSize);
            double[] imag = Arrays.copyOfRange(y, indx, indx + windowSize);
            for(int k = 0; k < real.length; k++) {
                System.out.print(real[k] + " ");
            }
            fft.fft(real, imag);

            System.out.println("");

            for(int k = 0; k < real.length; k++) {
                System.out.print(real[k] + " ");
            }

            currentWindow += 1;
            indx += windowSize;
            realRes.add(real);
            imagRes.add(imag);
        }

        ArrayList<ArrayList<double[]>> results = new ArrayList<>();
        results.add(realRes);
        results.add(imagRes);
        return results;
    }

    // Test
    public static void main(String[] args) {
        STFT stft = new STFT(16);

        int numWindows = 16;

        int N = (int) Math.pow(2,24);
        N = 4096;
        //N = 16384;
        //N = 262144;
        //N = 4194304;
        //N = N * 2;


        double[] re = new double[N];
        double[] im = new double[N];

        /*
        //Impulse stream
        for(int i = 0; i < N; i++) {
            if(i % 256 == 0) {
                re[i] = 2;
                im[i] = 0;
            } else if (i % 128 == 0) {
                re[i] = 1;
                im[i] = 0;
            } else
                re[i] = im[i] = 0;
        } */

        // Sin
        for(int i = 0; i < N; i++) {
            re[i] = Math.cos(2*Math.PI*i / (256 / 8));
            im[i] = 0;
        }

        long startTime = System.nanoTime();
        ArrayList<ArrayList<double[]>> freq =  stft.stft(re, im);
        long endTime = System.nanoTime();
        long execTime = endTime - startTime;

        System.out.println(execTime / (1000000) + " milliseconds");


        int len = freq.get(0).size();
        System.out.println("Real: ");
        for(int i = 0; i < len; i++) {
            System.out.print("[ ");
            for(int k = 0; k < freq.get(0).get(i).length; k++) {
                System.out.print(freq.get(0).get(i)[k] + " ");
            }
            System.out.println("]");
        }

        System.out.println("Imag: ");
        for(int i = 0; i < len; i++) {
            System.out.print("[ ");
            for(int k = 0; k < freq.get(1).get(i).length; k++) {
                System.out.print(freq.get(1).get(i)[k] + " ");
            }
            System.out.println("]");
        }
    }
}
