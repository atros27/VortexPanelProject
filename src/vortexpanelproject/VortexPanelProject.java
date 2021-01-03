/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vortexpanelproject;
//import smile.plot.Contour;
//import smile.plot.Heatmap;
//import smile.plot.PlotCanvas;

import javax.swing.*;
import java.awt.*;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.geom.Line2D;
//import java.io.File;
//import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

public class VortexPanelProject extends JPanel {

    public static Polygon airfoil;
    public static double[] airfoil_x_real, airfoil_y_real;
    public static double[] gamma;
    public static ArrayList<Line2D> streamlines;

    public static void main(String[] args) {
        // Prompt user input
        JTextField nacafield = new JTextField();
        JTextField chordfield = new JTextField();
        JTextField alphafield = new JTextField();
        JTextField vfield = new JTextField();
        JTextField pointsfield = new JTextField();
        Object[] fields = {
            "NACA Code", nacafield,
            "Chord [m]", chordfield,
            "Angle of Attack [deg]", alphafield,
            "Windspeed [m/s]", vfield,
            "# of points", pointsfield
        };
        JOptionPane.showConfirmDialog(null, fields, "Vortex Panel Project", JOptionPane.OK_CANCEL_OPTION);

        // Create airfoil polygon
        int naca = Integer.parseInt(nacafield.getText());
        double chord = Double.parseDouble(chordfield.getText());
        double AoA = Double.parseDouble(alphafield.getText());
        double v_inf = Double.parseDouble(vfield.getText());
        int num_points = Integer.parseInt(pointsfield.getText());
        int canvas_width = 900; //To be to-scale, must be 3x2
        int canvas_height = 600;
        double resolution = .01 * chord; //Field resolution in m
        airfoil = createAirfoil(naca, chord, AoA, num_points, canvas_width, canvas_height);
        gamma = panelMethod(v_inf, AoA);
        double[][] psi = createField(v_inf, chord, resolution);
        double[] chosen_lines = {-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5};
        streamlines = contour(psi, chosen_lines, canvas_width, canvas_height, chord, resolution);

        double[][] matrix = {{1, 4, 7}, {0, 3, 6}, {9, 2, 5}};
        double[] vector = {5, 10, 15};
        double[] solution = linearAlgebra(matrix, vector);
        for (int i = 0; i < solution.length; i++) {
            System.out.println(Double.toString(solution[i]));
        }

        JFrame frame = new JFrame();
        frame.setTitle("Vortex Panel Project");
        frame.setSize(canvas_width, canvas_height);
        Container contentPane = frame.getContentPane();
        contentPane.add(new VortexPanelProject());
        frame.setVisible(true);

        frame.addWindowListener(
                new WindowAdapter() {
            public void windowClosing(WindowEvent e) {
                System.exit(0);
            }
        });
    }

    static double[] backwardSubstitution(double[][] U, double[] y) {
        double[] x = new double[y.length];
        int i, j;
        double alpha;
        for (i = y.length - 1; i >= 0; i--) {
            alpha = 0;
            for (j = i + 1; j < y.length; j++) {
                alpha += U[i][j] * x[j];
            }
            x[i] = (y[i] - alpha) / U[i][i];
        }
        return x;
    }

    static ArrayList<Line2D> contour(double[][] z, double[] isolines, int canvas_width, int canvas_height, double chord, double resolution) {
        boolean[][] boolmat = new boolean[z.length][z[0].length];
        int a, b, c, d, id;
        ArrayList<Line2D> lines = new ArrayList();
        //create x-y coordinate matrix for ij referencing
        double[] x = new double[z[0].length];
        double[] y = new double[z.length];
        for (int j = 0; j < x.length; j++) {
            if (j == 0) {
                x[j] = -1 * chord;
            } else {
                x[j] = x[j - 1] + resolution;
            }
        }
        for (int i = 0; i < y.length; i++) {
            if (i == 0) {
                y[i] = chord;
            } else {
                y[i] = y[i - 1] - resolution;
            }
        }
        for (double iso : isolines) {
            for (int i = 0; i < z.length; i++) {
                for (int j = 0; j < z[0].length; j++) {
                    boolmat[i][j] = z[i][j] > iso;
                }
            }
            for (int i = 0; i < z.length - 1; i++) {
                for (int j = 0; j < z[0].length - 1; j++) {
                    a = boolmat[i][j] ? 1 : 0;
                    b = boolmat[i][j + 1] ? 1 : 0;
                    c = boolmat[i + 1][j + 1] ? 1 : 0;
                    d = boolmat[i + 1][j] ? 1 : 0;
                    id = a * 8 + b * 4 + c * 2 + d;

                    double dx1, dx2, dy1, dy2;
                    Line2D line1, line2;

                    switch (id) {
                        case 1:
                        case 14:
                            //a-d vert interp
                            dy1 = interpolate(z[i][j] - iso, z[i + 1][j] - iso, resolution);
                            //c-d horz interp
                            dx2 = interpolate(z[i + 1][j] - iso, z[i + 1][j + 1] - iso, resolution);
                            //create line (and add to list)
                            line1 = pointsToLine(x[j], y[i] - dy1, x[j] + dx2, y[i + 1], chord, canvas_width, canvas_height);
                            lines.add(line1);
                            break;
                        case 2:
                        case 13:
                            //c-d horz interp
                            dx1 = interpolate(z[i + 1][j] - iso, z[i + 1][j + 1] - iso, resolution);
                            //b-c vert interp
                            dy2 = interpolate(z[i][j + 1] - iso, z[i + 1][j + 1] - iso, resolution);
                            //create line
                            line1 = pointsToLine(x[j] + dx1, y[i], x[j + 1], y[i] - dy2, chord, canvas_width, canvas_height);
                            lines.add(line1);
                            break;
                        case 3:
                        case 12:
                            //a-d vert interp
                            dy1 = interpolate(z[i][j] - iso, z[i + 1][j] - iso, resolution);
                            //b-c vert interp
                            dy2 = interpolate(z[i][j + 1] - iso, z[i + 1][j + 1] - iso, resolution);
                            //create line
                            line1 = pointsToLine(x[j], y[i] - dy1, x[j + 1], y[i] - dy2, chord, canvas_width, canvas_height);
                            lines.add(line1);
                            break;
                        case 4:
                        case 11:
                            //a-b horz interp
                            dx1 = interpolate(z[i][j] - iso, z[i][j + 1] - iso, resolution);
                            //b-c vert interp
                            dy2 = interpolate(z[i][j + 1] - iso, z[i + 1][j + 1] - iso, resolution);
                            //create line
                            line1 = pointsToLine(x[j] + dx1, y[i], x[j + 1], y[i] - dy2, chord, canvas_width, canvas_height);
                            lines.add(line1);
                            break;
                        case 5:
                            //a-d vert interp
                            dy1 = interpolate(z[i][j] - iso, z[i + 1][j] - iso, resolution);
                            //a-b horz interp
                            dx2 = interpolate(z[i][j] - iso, z[i][j + 1] - iso, resolution);
                            //create line
                            line1 = pointsToLine(x[j], y[i] - dy1, x[j] + dx2, y[i], chord, canvas_width, canvas_height);
                            lines.add(line1);

                            //b-c vert interp
                            dy1 = interpolate(z[i][j + 1] - iso, z[i + 1][j + 1] - iso, resolution);
                            //c-d horz interp
                            dx2 = interpolate(z[i + 1][j] - iso, z[i + 1][j + 1] - iso, resolution);
                            //create line
                            line2 = pointsToLine(x[j] + dx2, y[i], x[j + 1], y[i] - dy1, chord, canvas_width, canvas_height);
                            lines.add(line2);
                            break;
                        case 6:
                        case 9:
                            //a-b horz interp
                            dx1 = interpolate(z[i][j] - iso, z[i][j + 1] - iso, resolution);
                            //c-d horz interp
                            dx2 = interpolate(z[i + 1][j] - iso, z[i + 1][j + 1] - iso, resolution);
                            //create line
                            line1 = pointsToLine(x[j] + dx1, y[i], x[j] + dx2, y[i + 1], chord, canvas_width, canvas_height);
                            lines.add(line1);
                            break;
                        case 7:
                        case 8:
                            //a-b horz interp
                            dx1 = interpolate(z[i][j] - iso, z[i][j + 1] - iso, resolution);
                            //a-d vert interp
                            dy2 = interpolate(z[i][j] - iso, z[i + 1][j] - iso, resolution);
                            //create line
                            line1 = pointsToLine(x[j] + dx1, y[i], x[j], y[i] - dy2, chord, canvas_width, canvas_height);
                            lines.add(line1);
                            break;
                        case 10:
                            //a-b horz interp
                            dx1 = interpolate(z[i][j] - iso, z[i][j + 1] - iso, resolution);
                            //b-c vert interp
                            dy2 = interpolate(z[i][j + 1] - iso, z[i + 1][j + 1] - iso, resolution);
                            //create line
                            line1 = pointsToLine(x[j] + dx1, y[i], x[j + 1], y[i] - dy2, chord, canvas_width, canvas_height);
                            lines.add(line1);

                            //c-d horz interp
                            dx1 = interpolate(z[i + 1][j] - iso, z[i + 1][j + 1] - iso, resolution);
                            //a-d vert interp
                            dy2 = interpolate(z[i][j] - iso, z[i + 1][j] - iso, resolution);
                            //create line
                            line2 = pointsToLine(x[j] + dx1, y[i + 1], x[j], y[i] - dy2, chord, canvas_width, canvas_height);
                            lines.add(line2);
                            break;
                        default: // No contour line
                            break;
                    }
                }
            }
        }
        return lines;
    }

    static Polygon createAirfoil(int naca, double chord, double AoA, int num_points, int width, int height) {
        //break NACA into parts
        double m = (double) (naca / 1000) / 100;
        double p = (double) ((naca % 1000) / 100) / 10;
        double t = (double) (naca % 100) / 100;

        //define initial x
        double[] x = new double[num_points];
        int i;
        double step_size = chord / (num_points / 2);
        x[0] = chord;
        for (i = 1; i < num_points; i++) {
            if (i < num_points / 2 + 1) {
                x[i] = x[i - 1] - step_size;
            } else {
                x[i] = x[i - 1] + step_size;
            }
        }

        //define final x and y
        double yt;
        double yc;
        double theta;
        double dyc;//intermediate vars for theta
        double dx;
        double[] y = new double[num_points];
        for (i = 0; i < num_points; i++) {
            //define intermediate variables
            yt = 5 * t * chord * (.2969 * Math.sqrt(Math.abs(x[i] / chord)) - .126 * (x[i] / chord) - .3516 * Math.pow(x[i] / chord, 2) + .2843 * Math.pow(x[i] / chord, 3) - .1015 * Math.pow(x[i] / chord, 4));
            if (x[i] < p * chord) {
                yc = m / Math.pow(p, 2) * (2 * p * (x[i] / chord) - Math.pow(x[i] / chord, 2));
                dyc = 2 * m * (p - x[i] / chord);
                dx = Math.pow(p, 2);
                theta = Math.atan2(dyc, dx);
            } else {
                yc = m / Math.pow(1 - p, 2) * (1 - 2 * p + 2 * p * (x[i] / chord) - Math.pow(x[i] / chord, 2));
                dyc = 2 * m * (p - x[i] / chord);
                dx = Math.pow(1 - p, 2);
                theta = Math.atan2(dyc, dx);
            }

            //create finalized x and y
            if (i < num_points / 2 + 1) { //y >= 0
                x[i] -= yt * Math.sin(theta);
                y[i] = yc + yt * Math.cos(theta);
            } else { //y < 0
                x[i] += yt * Math.sin(theta);
                y[i] = yc - yt * Math.cos(theta);
            }
        }

        int[] xpixels = new int[num_points];
        int[] ypixels = new int[num_points];
        for (i = 0; i < num_points; i++) {
            //rotate points
            x[i] = x[i] * Math.cos(-1 * AoA / 180 * Math.PI) - y[i] * Math.sin(-1 * AoA / 180 * Math.PI);
            y[i] = x[i] * Math.sin(-1 * AoA / 180 * Math.PI) + y[i] * Math.cos(-1 * AoA / 180 * Math.PI);

            //transform points to pixels
            xpixels[i] = (int) ((width / 3) * (1 + x[i] / chord));
            ypixels[i] = (int) ((height / 2) * (1 - y[i] / chord));
        }
        airfoil_x_real = x;
        airfoil_y_real = y;

        //Create Polygon
        Polygon ans = new Polygon(xpixels, ypixels, num_points);
        return ans;
    }

    static double[][] createField(double v_inf, double chord, double resolution) {
        double[] x = new double[3 * ((int) (chord / resolution)) + 1];
        double[] y = new double[2 * ((int) (chord / resolution)) + 1];
        double[][] z = new double[y.length][x.length];
        for (int i = 0; i < y.length; i++) {
            if (i == 0) {
                y[i] = chord;
            } else {
                y[i] = y[i - 1] - resolution;
            }
        }
        for (int j = 0; j < x.length; j++) {
            if (j == 0) {
                x[j] = -1 * chord;
            } else {
                x[j] = x[j - 1] + resolution;
            }
        }

        double r;
        for (int i = 0; i < y.length; i++) {
            for (int j = 0; j < x.length; j++) {
                z[i][j] = v_inf * y[i];
                for (int k = 0; k < gamma.length; k++) {
                    r = Math.sqrt(Math.pow(x[j] - airfoil_x_real[k % airfoil_x_real.length], 2) + Math.pow(y[i] - airfoil_y_real[k % airfoil_y_real.length], 2));
                    z[i][j] += gamma[k] * Math.log(r);
                    //System.out.printf("r: %f, lnr: %f, gamma: %f\n", r, Math.log(r), gamma[k]);
                }
                //System.out.printf("%f ",z[i][j]);
            }
            //System.out.printf("\n");
        }
        return z;
    }

    static double[][] decompose(double[][] A) {
        double ratio;
        double[][] LU = Arrays.stream(A).map(e -> e.clone()).toArray($ -> A.clone());
        for (int i = 0; i < LU.length - 1; i++) {
            for (int j = i + 1; j < LU.length; j++) {
                ratio = LU[j][i] / LU[i][i];
                for (int k = i; k < LU.length; k++) {
                    LU[j][k] -= ratio * LU[i][k];
                }
                LU[j][i] = ratio;
            }
        }

        //for (int i = 0; i < LU.length; i++) {
        //    for (int j = 0; j < LU.length; j++) {
        //        System.out.printf("%f ", LU[i][j]);
        //    }
        //    System.out.printf("\n");
        //}

        return LU;
    }

    static double[] forwardSubstitution(double[][] L, double[] b) {
        double[] y = new double[b.length];
        double alpha;
        for (int i = 0; i < L.length; i++) {
            alpha = 0;
            for (int j = 0; j < i; j++) {
                alpha += L[i][j] * y[j];
            }
            y[i] = b[i] - alpha;
        }
        
        for(int i = 0; i < L.length; i++) {
            System.out.printf("%f\n", y[i]);
        }
        
        return y;
    }

    static double[] linearAlgebra(double[][] A, double[] b) {
        double[] x = new double[b.length];
        if (A.length != A[0].length) {
            System.out.println("Error: Matrix is not square");
        } else if (A[0].length != b.length) {
            System.out.println("Error: Matrix and vector incompatible");
        } else {
            double[][] LU = decompose(A);
            double[] y = forwardSubstitution(LU, b);
            x = backwardSubstitution(LU, y);
        }
        return x;
    }

    static double interpolate(double y1, double y2, double resolution) {  // find portion along resolution (from point 1) that crosses isoline
        double m = (y2 - y1) / resolution;
        double x = -1 * y1 / m;
        return x;
    }

    @Override
    public void paintComponent(Graphics g) {
        super.paintComponent(g); // Paint Airfoil
        g.drawPolygon(airfoil);
        g.setColor(Color.BLACK);
        g.fillPolygon(airfoil);
        Graphics2D g2 = (Graphics2D) g;
        g2.setColor(Color.LIGHT_GRAY);
        for (int i = 0; i < streamlines.size(); i++) {
            g2.draw(streamlines.get(i));
        }
    }

    static double[] panelMethod(double v_inf, double AoA) {
        AoA *= Math.PI / 180;
        double[][] AN = new double[airfoil_x_real.length + 1][airfoil_x_real.length + 1];
        double[] RHS = new double[airfoil_x_real.length + 1];
        double[] theta = new double[airfoil_x_real.length];
        double[] S = new double[airfoil_x_real.length];
        for (int i = 0; i < airfoil_x_real.length; i++) {
            theta[i] = Math.atan2(airfoil_y_real[(i + 1) % airfoil_y_real.length] - airfoil_y_real[i],
                    airfoil_x_real[(i + 1) % airfoil_x_real.length] - airfoil_x_real[i]);
            S[i] = Math.sqrt(Math.pow(airfoil_y_real[(i + 1) % airfoil_y_real.length] - airfoil_y_real[i], 2)
                    + Math.pow(airfoil_x_real[(i + 1) % airfoil_x_real.length] - airfoil_x_real[i], 2));
        }
        double A, B, C, D, E, F, G, P, Q, Cn1, Cn2, Cn2_prev;
        for (int i = 0; i < airfoil_x_real.length + 1; i++) {
            if (i == airfoil_x_real.length) {
                RHS[i] = 0;
            } else {
                RHS[i] = Math.sin(theta[i] - AoA);
            }
            Cn2_prev = 0;
            for (int j = 0; j < airfoil_x_real.length + 1; j++) { // AN Matrix
                if (i == airfoil_x_real.length) { // Last row of AN
                    if (j == 0 || j == airfoil_x_real.length) {
                        AN[i][j] = 1;
                    } else {
                        AN[i][j] = 0;
                    }
                } else if (i == j) { // Diagonal
                    Cn1 = -1;
                    Cn2 = 1;
                    AN[i][j] = Cn1 + Cn2_prev;
                    Cn2_prev = Cn2;
                } else { // Non-Diagonal
                    if (j == airfoil_x_real.length) {
                        AN[i][j] = Cn2_prev;
                    } else {
                        A = -1 * ((airfoil_x_real[i] - airfoil_x_real[j]) * Math.cos(theta[j]) + (airfoil_y_real[i] - airfoil_y_real[j]) * Math.sin(theta[j]));
                        B = Math.pow(airfoil_x_real[i] - airfoil_x_real[j], 2) + Math.pow(airfoil_y_real[i] - airfoil_y_real[j], 2);
                        C = Math.sin(theta[i] - theta[j]);
                        D = Math.cos(theta[i] - theta[j]);
                        E = (airfoil_x_real[i] - airfoil_x_real[j]) * Math.sin(theta[j]) - (airfoil_y_real[i] - airfoil_y_real[j]) * Math.cos(theta[j]);
                        F = Math.log(1 + (Math.pow(S[j], 2) + 2 * A * S[j]) / B);
                        G = Math.atan2(E * S[j], B + A * S[j]);
                        P = (airfoil_x_real[i] - airfoil_x_real[j]) * Math.sin(theta[i] - 2 * theta[j]) + (airfoil_y_real[i] - airfoil_y_real[j]) * Math.cos(theta[i] - 2 * theta[j]);
                        Q = (airfoil_x_real[i] - airfoil_x_real[j]) * Math.cos(theta[i] - 2 * theta[j]) - (airfoil_y_real[i] - airfoil_y_real[j]) * Math.sin(theta[i] - 2 * theta[j]);
                        Cn2 = D + .5 * Q * F / S[j] - (A * C + D * E) * G / S[j];
                        Cn1 = .5 * D * F + C * G - Cn2;
                        AN[i][j] = Cn1 + Cn2_prev;
                        Cn2_prev = Cn2;
                    }
                }
            }
        }
        
        for(int i = 0; i<AN.length; i++) {
            for(int j = 0; j<AN.length; j++) {
                System.out.printf("%f ", AN[i][j]);
            }
            System.out.printf("\n");
        }
        
        double[] gamma = linearAlgebra(AN, RHS);
        for (int i = 0; i < gamma.length; i++) {
            gamma[i] *= v_inf;
        }
        return gamma;
    }

    static Line2D pointsToLine(double x1, double y1, double x2, double y2, double chord, int width, int height) {
        int px1 = (int) ((width / 3) * (1 + x1 / chord));
        int px2 = (int) ((width / 3) * (1 + x2 / chord));
        int py1 = (int) ((height / 2) * (1 - y1 / chord));
        int py2 = (int) ((height / 2) * (1 - y2 / chord));
        return new Line2D.Double(px1, px2, py1, py2);
    }
}
