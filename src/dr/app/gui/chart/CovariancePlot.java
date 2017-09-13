/*
 * CovariancePlot.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.gui.chart;

import dr.stats.DiscreteStatistics;
import dr.stats.Variate;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;

/**
 * Description:	A covariance ellipse plot.
 *
 * @author Guy Baele
 */
public class CovariancePlot extends Plot.AbstractPlot {

    private final int POSITIVE_CORRELATION_DEGREE = 35;
    private final int NEGATIVE_CORRELATION_DEGREE = 235;

    //colors from ColorBrewer2.org
    private final Color STRONG_POSITIVE_CORRELATION = new Color(227,74,51);
    private final Color MODERATE_POSITIVE_CORRELATION = new Color(253,187,132);
    private final Color WEAK_POSITIVE_CORRELATION = new Color(254,232,200);

    private final Color STRONG_NEGATIVE_CORRELATION = new Color(49,130,189);
    private final Color MODERATE_NEGATIVE_CORRELATION = new Color(158,202,225);
    private final Color WEAK_NEGATIVE_CORRELATION = new Color(222,235,247);

    private int plotCount;

    public CovariancePlot(java.util.List<Double> data, int minimumBinCount) {
        super(data, data);
        System.out.println("CovariancePlot: " + data.size());
        setName("null");
        this.plotCount = 1;
    }

    public CovariancePlot(java.util.List<Double> xData, java.util.List<Double> yData) {
        super(xData, yData);
        System.out.println("xData: " + xData.size());
        System.out.println("yData: " + yData.size());
        setName("null");
        this.plotCount = 1;
    }

    public CovariancePlot(String name, java.util.List<Double> xData, java.util.List<Double> yData) {
        super(xData, yData);
        setName(name);
        this.plotCount = 1;
    }

    public void setTotalPlotCount(int count) {
        this.plotCount = count;
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {

        //System.out.println("CovariancePlot: paintData");
        //System.out.println("PlotNumber = " + plotNumber);
        //System.out.println("TotalPlotCount = " + this.plotCount);
        int xCount = xData.getCount();
        int yCount = yData.getCount();
        double[] xDataArray = new double[xCount];
        double[] yDataArray = new double[yCount];
        for (int i = 0; i < xCount; i++) {
            xDataArray[i] = (Double)xData.get(i);
            yDataArray[i] = (Double)yData.get(i);
        }

        double covariance = DiscreteStatistics.covariance(xDataArray, yDataArray);
        double stdevX = DiscreteStatistics.stdev(xDataArray);
        double stdevY = DiscreteStatistics.stdev(yDataArray);

        double correlation = covariance / (stdevX * stdevY);

        System.out.println("plotNumber: " + plotNumber + " ; correlation = " + correlation);

        Color fillColor;
        if (correlation >= 0.0 && correlation < 0.30) {
            fillColor = WEAK_POSITIVE_CORRELATION;
        } else if (correlation >= 0.30 && correlation < 0.70) {
            fillColor = MODERATE_POSITIVE_CORRELATION;
        } else if (correlation >= 0.70 && correlation <= 1.0) {
            fillColor = STRONG_POSITIVE_CORRELATION;
        } else if (correlation < 0.0 && correlation >= -0.30) {
            fillColor = WEAK_NEGATIVE_CORRELATION;
        } else if (correlation < -0.30 && correlation >= -0.70) {
            fillColor = MODERATE_NEGATIVE_CORRELATION;
        } else {
            fillColor = STRONG_NEGATIVE_CORRELATION;
        }
        g2.setColor(fillColor);

        //g2.setPaint(linePaint);
        //g2.setStroke(lineStroke);

        double x1 = (plotNumber/(int)Math.sqrt(plotCount))+0.65;
        double y1 = (plotNumber%(int)Math.sqrt(plotCount))+0.65;
        double x2 = (plotNumber/(int)Math.sqrt(plotCount))+1.35;
        double y2 = (plotNumber%(int)Math.sqrt(plotCount))+1.35;

        //System.out.println("(" + x1 + "," + y1 + ") (" + x2 + "," + y2 + ")");

        drawRect(g2, x1, y1, x2, y2);
        g2.drawString("plot" + plotNumber, (float)transformX(x1), (float)transformY(y1));

        double rotationDegree = NEGATIVE_CORRELATION_DEGREE;
        if (correlation > 0) {
            rotationDegree = POSITIVE_CORRELATION_DEGREE;
        }

        //TODO alter ellipse height according to correlation
        double ellipseWidth = Math.abs(transformX(x2)-transformX(x1));
        double ellipseHeight = 0.5*Math.abs(transformY(y2) - transformY(y1));

        AffineTransform oldTransform = g2.getTransform();

        Shape ellipse = new Ellipse2D.Double(transformX(x1), transformY(y2), ellipseWidth, ellipseHeight);
        //this transformation rotates around the upper right rectangle corner
        //Shape rotatedEllipse = AffineTransform.getRotateInstance(rotationDegree, transformX(x1), transformY(y2)).createTransformedShape(ellipse);

        //rotate around rectangle center
        Shape rotatedEllipse = AffineTransform.getRotateInstance(rotationDegree, transformX((x1+x2)/2.0), transformY((y1+y2)/2.0)).createTransformedShape(ellipse);
        //rotate around ellipse center
        //Shape rotatedEllipse = AffineTransform.getRotateInstance(rotationDegree, transformX((x1+x2)/2.0), transformY((y1+y2)/2.0)).createTransformedShape(ellipse);

        g2.fill(rotatedEllipse);

        //draw rotation center in black
        g2.setColor(Color.BLACK);
        g2.fill(new Ellipse2D.Double(transformX((x1+x2)/2.0), transformY((y1+y2)/2.0), 5, 5));

        //draw potential new rotation center
        //g2.fill(new Ellipse2D.Double(transformX((x1+x2)/2.0), transformY((y1+y2)/2.0), 5, 5));

        g2.setTransform(oldTransform);

    }

}
