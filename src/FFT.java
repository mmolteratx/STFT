import java.lang.Math.*;

public class FFT {

    int n, m;

    double[] cos;
    double[] sin;

    double[] window;

    public FFT(int n) {
        this.n = n;
        this.m = (int)(Math.log(n) / Math.log(2));

        if(n != (1<<m)) {
            throw new RuntimeException("FFT length must be power of 2");
        }

        cos = new double[n/2];
        sin = new double[n/2];

        for(int i = 0; i < n/2; i++) {
            cos[i] = Math.cos(-2 * Math.PI * i/n);
            sin[i] = Math.sin(-2 * Math.PI * i/n);
        }

        makeWindow();
    }

    // TODO: Try different windows
    protected void makeWindow() {
        // Blackman window
        window = new double[n];
        for(int i = 0; i < window.length; i++) {
            window[i] = 0.42 - 0.5 * Math.cos(2*Math.PI*i/(n-1)) +
                    0.08 * Math.cos(4*Math.PI*i/(n-1));
        }
    }

    public double[] getWindow() {
        return window;
    }

    public void fft(double[] x, double[] y) {

        int i, j, k, n1, n2, a;
        double c, s, e, t1, t2;

        // Bit reverse
        j = 0;
        n2 = n / 2;
        for(i = 1; i < n - 1; i++) {
            n1 = n2;
            while(j >= n1) {
                j = j - n1;
                n1 = n1/2;
            }
            j = j + n1;

            if(i < j) {
                t1 = x[i];
                x[i] = x[j];
                t1 = y[i];
                y[i] = y[j];
                y[j] = t1;
            }
        }

        // FFT
        n1 = 0;
        n2 = 1;
        for(i = 0; i < m; i++) {
            n1 = n2;
            n2 = n2 + n2;
            a = 0;

            for(j = 0; j < n1; j++) {
                c = cos[a];
                s = sin[a];
                a += 1 << (m-i-1);

                for(k = j; k < n; k = k + n2) {
                    t1 = c*x[k + n1] - s*y[k + n1];
                    t2 = s*x[k + n1] - c*y[k + n1];
                    x[k + n1]  = x[k] - t1;
                    y[k + n1] = y[k] - t2;
                    x[k] = x[k] + t1;
                    y[k] = y[k] + t2;
                }
            }
        }
    }

    public static void main(String[] args) {

        int N = (int) Math.pow(2,24);

        FFT fft = new FFT(N);

        double[] re = new double[N];
        double[] im = new double[N];

        /*
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
        } */

        // Sin
        for(int i = 0; i < N; i++) {
            re[i] = Math.cos(2*Math.PI*i / (N / 8));
            im[i] = 0;
        }
        //beforeAfter(fft, re, im);

        long startTime = System.nanoTime();
        fft.fft(re, im);
        long endTime = System.nanoTime();
        long execTime = endTime - startTime;

        System.out.println(execTime / (1000000) + " milliseconds");
    }

    protected static void beforeAfter(FFT fft, double[] re, double[] im) {
        System.out.println("Before: ");
        printReIm(re, im);
        fft.fft(re, im);
        System.out.println("After: ");
        printReIm(re, im);
    }

    protected static void printReIm(double[] re, double[] im) {
        System.out.print("Re: [");
        for (int i = 0; i < re.length; i++) {
            System.out.print(((int) (re[i] * 1000) / 1000.0) + " ");
        }
        System.out.print("]\nIm: [");
        for (int i = 0; i < im.length; i++) {
            System.out.print(((int) (im[i] * 1000) / 1000.0) + " ");
        }
        System.out.println("]");
    }
}
