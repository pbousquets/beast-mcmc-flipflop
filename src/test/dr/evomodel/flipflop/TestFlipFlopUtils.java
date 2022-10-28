package test.dr.evomodel.flipflop;
import java.util.Arrays;
import dr.util.FlipFlopUtils;

import junit.framework.TestCase;

/**
 * Test the flipflop utils
 *
 * @author Pablo Bousquets
 */
public class TestFlipFlopUtils extends TestCase {

    interface Instance {

        double[] getSequence();

        int getStart();

        int getEnd();

        int getNum();

        double[] getLinspaceResult();

        int[] getDigitizedResult();

        double[] getRemapResult();
    }

    Instance test = new Instance() {
        public double[] getSequence(){
            return new double[]{
                    0.1, 0.4, 0.14, 0.12, 0.99, .05, 0.01
            };
        };

        public int getStart() {
            return 0;
        }

        public int getEnd() {
            return 1;
        }

        public int getNum() {
            return 10;
        }

        public double[] getLinspaceResult() {
            return new double[]{
                    0.        , 0.11111111, 0.22222222, 0.33333333, 0.44444444,
                    0.55555556, 0.66666667, 0.77777778, 0.88888889, 1
            };
        }
        public int[] getDigitizedResult() {
            return new int[]{
                    1, 4, 2, 2, 9, 1, 1
            };
        }

        public double[] getRemapResult() {
            return new double[]{
                    0.1, 0.4, 0.2, 0.2, 0.9, 0.1, 0.1
            };
        }
    };

    public void testLinspace(){
        double[] bins = FlipFlopUtils.linspace(test.getStart(), test.getEnd(), test.getNum());
        System.out.println(Arrays.toString(bins));
        System.out.println(Arrays.toString(test.getLinspaceResult()));
        assertEquals(bins, test.getLinspaceResult(), 1e-4);
    }

    public void testDigitize(){
        double[] bins = FlipFlopUtils.linspace(test.getStart(), test.getEnd(), test.getNum());
        int[] digitized = FlipFlopUtils.digitize(test.getSequence(), bins);
        System.out.println(Arrays.toString(digitized));
        System.out.println(Arrays.toString(test.getDigitizedResult()));
        assertEquals(digitized, test.getDigitizedResult(), 1e-15);

    }

    public void testMapBack(){
        double[] bins = FlipFlopUtils.linspace(test.getStart(), test.getEnd(), test.getNum());
        int[] digitized = FlipFlopUtils.digitize(test.getSequence(), bins);
        double[] remapped = FlipFlopUtils.remapDigitized(digitized, test.getNum());
        System.out.println(Arrays.toString(remapped));
        System.out.println(Arrays.toString(test.getRemapResult()));
        assertEquals(remapped, test.getRemapResult(), 1e-4);
    }

    protected void assertEquals(double[] a, double[] b, double accuracy) {
        assertEquals(a.length, b.length);
        for (int i = 0; i < a.length; i++) {
            assertEquals(a[i], b[i], accuracy);
        }
    }

    protected void assertEquals(int[] a, int[] b, double accuracy) {
        assertEquals(a.length, b.length);
        for (int i = 0; i < a.length; i++) {
            assertEquals(a[i], b[i], accuracy);
        }
    }
}
