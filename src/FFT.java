import java.lang.Math.*;

public class FFT {

    int n, m;

    double[] cos;
    double[] sin;

    public FFT(int n) {
        this.n = n;
        this.m = (int)(Math.log(n) / Math.log(2));

        if(n != (1 << m)) {
            throw new RuntimeException("FFT length must be power of 2");
        }

        cos = new double[n >> 1];
        sin = new double[n >> 1];

        for(int i = 0; i < n/2; i++) {
            cos[i] = Math.cos(-2 * Math.PI * i/n);
            sin[i] = Math.sin(-2 * Math.PI * i/n);
        }

    }

    public void fft(double[] x, double[] y) {

        int i, j, k, n1, n2, a;
        double c, s, temp1, temp2;

        // Bit reverse
        j = 0;
        n2 = n >> 1;
        for(i = 1; i < n - 1; i++) {
            n1 = n2;
            while(j >= n1) {
                j -= n1;
                n1 = n1 >> 1;
            }
            j += n1;

            if(i < j) {
                temp1 = x[i];
                x[i] = x[j];
                x[j] = temp1;
                temp1 = y[i];
                y[i] = y[j];
                y[j] = temp1;
            }
        }

        // FFT
        n2 = 1;
        // O(logn)
        for(i = 0; i < m; i++) {
            n1 = n2;
            n2 = n2 << 1;
            a = 0;

            // O(n)
            for(j = 0; j < n1; j++) {
                c = cos[a];
                s = sin[a];
                a += 1 << (m-i-1);

                 
                for(k = j; k < n; k += n2) {
                    temp1 = c * x[k+n1] - s * y[k+n1];
                    temp2 = s * x[k+n1] + c * y[k+n1];
                    x[k+n1] = x[k] - temp1;
                    y[k+n1] = y[k] - temp2;
                    x[k] += temp1;
                    y[k] += temp2;
                }
            }
        }
    }

    public static void main(String[] args) {

        int N = (int) Math.pow(2,24);
        N = 256;
        //N = 16384;
        //N = 262144;
        //N = 4194304;
        //N = N * 2;

        FFT fft = new FFT(N);

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
            re[i] = Math.cos(2*Math.PI*i / (N / 8));
            im[i] = 0;
        }

        long startTime = System.nanoTime();
        fft.fft(re, im);
        long endTime = System.nanoTime();
        long execTime = endTime - startTime;

        System.out.println(execTime / (1000000) + " milliseconds");

        System.out.print("Real: [");
        for (int i = 0; i < re.length; i++) {
            System.out.print(((int) (re[i] * 1000) / 1000.0) + ",");
        }
        System.out.print("]\nImaginary: [");
        for (int i = 0; i < im.length; i++) {
            System.out.print(((int) (im[i] * 1000) / 1000.0) + " ");
        }
        System.out.println("]");

    }
}
