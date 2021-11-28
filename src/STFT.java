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

        while(currentWindow < numWindows) {
            double[] real = Arrays.copyOfRange(x, indx, indx + windowSize);
            double[] imag = Arrays.copyOfRange(y, indx, indx + windowSize);
            FFT fft = new FFT(windowSize);
            fft.fft(real, imag);
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
        STFT stft = new STFT(3);

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

        ArrayList<ArrayList<double[]>> freq =  stft.stft(re, im);

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
