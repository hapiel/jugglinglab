// Animator.java
//
// Copyright 2018 by Jack Boyce (jboyce@gmail.com) and others

/*
    This file is part of Juggling Lab.

    Juggling Lab is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Juggling Lab is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Juggling Lab; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

package jugglinglab.core;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ResourceBundle;

import javax.imageio.*;
import javax.imageio.metadata.*;
import javax.imageio.stream.ImageOutputStream;
import javax.imageio.stream.MemoryCacheImageOutputStream;
import org.w3c.dom.Node;

import jugglinglab.jml.*;
import jugglinglab.renderer.Renderer2D;
import jugglinglab.util.*;

// import gifwriter.GIFAnimWriter;


public class Animator {

    protected JMLPattern		pat;
    protected AnimationPrefs	jc;
    protected jugglinglab.renderer.Renderer	ren1 = null, ren2 = null;
    protected Coordinate		overallmax = null, overallmin = null;

    protected int[]				animpropnum = null, temppropnum = null;
    protected Permutation		invpathperm = null;
    protected int				num_frames;
    protected double			sim_time;
    protected double			sim_interval_secs;
    protected long				real_interval_millis;

    protected double[]          camangle;
    protected double[]			camangle1;     // for stereo display
    protected double[]			camangle2;

    protected Dimension         dim;


    public Animator() {
        this.camangle = new double[2];
        this.camangle1 = new double[2];
        this.camangle2 = new double[2];
		jc = new AnimationPrefs();
    }


    public void restartAnimator(JMLPattern pat, AnimationPrefs newjc)
                    throws JuggleExceptionUser, JuggleExceptionInternal {
        // try to lay out new pattern first so that if there's an error we
        // won't stop the current animation
        if ((pat != null) && !pat.isLaidout())
            pat.layoutPattern();

        if (pat != null)	this.pat = pat;
        if (newjc != null)	this.jc = newjc;

        if (this.pat == null)
            return;

		ren1 = new Renderer2D();
		if (this.jc.stereo)
			ren2 = new Renderer2D();

        ren1.setPattern(this.pat);
        if (this.jc.stereo)
            ren2.setPattern(this.pat);

        initAnimator();

        if (this.pat.getNumberOfJugglers() == 1) {
            this.camangle[0] = JLMath.toRad(0.0);
            this.camangle[1] = JLMath.toRad(90.0);
        } else {
            this.camangle[0] = JLMath.toRad(340.0);
            this.camangle[1] = JLMath.toRad(70.0);
        }

        setCameraAngle(camangle);

        if (jugglinglab.core.Constants.DEBUG_LAYOUT)
            System.out.println(this.pat);
    }

    public Dimension getDimension() {
        return this.dim;
    }

    public void setDimension(Dimension dim) {
        this.dim = dim;
        if (ren1 != null)
            syncRenderersToSize();
    }

    public double[] getCameraAngle() {
        return this.camangle;
    }

    protected void setCameraAngle(double[] ca) {
        while (ca[0] < 0.0)
            ca[0] += JLMath.toRad(360.0);
        while (ca[0] >= JLMath.toRad(360.0))
            ca[0] -= JLMath.toRad(360.0);

        this.camangle[0] = ca[0];
        this.camangle[1] = ca[1];

        if (jc.stereo) {
            this.camangle1[0] = ca[0] - 0.05;
            this.camangle1[1] = ca[1];
            this.ren1.setCameraAngle(this.camangle1);
            this.camangle2[0] = ca[0] + 0.05;
            this.camangle2[1] = ca[1];
            ren2.setCameraAngle(this.camangle2);
        } else {
            this.camangle1[0] = ca[0];
            this.camangle1[1] = ca[1];
            ren1.setCameraAngle(this.camangle1);
        }
    }



    public void drawFrame(double sim_time, Graphics g, boolean draw_axes) throws JuggleExceptionInternal {
        int[] pnum = this.animpropnum;

        if (this.jc.stereo) {
            this.ren1.drawFrame(sim_time, pnum,
                                g.create(0, 0, this.dim.width/2, this.dim.height));
            this.ren2.drawFrame(sim_time, pnum,
                                g.create(this.dim.width/2, 0, this.dim.width/2, this.dim.height));
        } else {
            this.ren1.drawFrame(sim_time, pnum, g);
        }

        if (draw_axes) {
            if (g instanceof Graphics2D) {
                Graphics2D g2 = (Graphics2D)g;
                g2.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
            }

            double[] ca = ren1.getCameraAngle();
            double theta = ca[0];
            double phi = ca[1];

            double xya = 30.0;
            double xyb = xya * Math.sin(90.0*0.0174532925194 - phi);
            double zlen = xya * Math.cos(90.0*0.0174532925194 - phi);
            int cx = 38;
            int cy = 45;
            int xx = cx + (int)(0.5 - xya * Math.cos(theta));
            int xy = cy + (int)(0.5 + xyb * Math.sin(theta));
            int yx = cx + (int)(0.5 - xya * Math.cos(theta +
                                                    90.0*0.0174532925194));
            int yy = cy + (int)(0.5 + xyb * Math.sin(theta +
                                                    90.0*0.0174532925194));
            int zx = cx;
            int zy = cy - (int)(0.5 + zlen);

            g.setColor(Color.green);
            g.drawLine(cx, cy, xx, xy);
            g.drawLine(cx, cy, yx, yy);
            g.drawLine(cx, cy, zx, zy);
            g.fillOval(xx-2, xy-2, 5, 5);
            g.fillOval(yx-2, yy-2, 5, 5);
            g.fillOval(zx-2, zy-2, 5, 5);
            g.drawString("x", xx-2, xy-4);
            g.drawString("y", yx-2, yy-4);
            g.drawString("z", zx-2, zy-4);
        }

        if (this.jc.stereo && draw_axes) {
            double[] ca = ren2.getCameraAngle();
            double theta = ca[0];
            double phi = ca[1];

            double xya = 30.0;
            double xyb = xya * Math.sin(90.0*0.0174532925194 - phi);
            double zlen = xya * Math.cos(90.0*0.0174532925194 - phi);
            int cx = 38 + this.dim.width/2;
            int cy = 45;
            int xx = cx + (int)(0.5 - xya * Math.cos(theta));
            int xy = cy + (int)(0.5 + xyb * Math.sin(theta));
            int yx = cx + (int)(0.5 - xya * Math.cos(theta + 90.0*0.0174532925194));
            int yy = cy + (int)(0.5 + xyb * Math.sin(theta + 90.0*0.0174532925194));
            int zx = cx;
            int zy = cy - (int)(0.5 + zlen);

            g.setColor(Color.green);
            g.drawLine(cx, cy, xx, xy);
            g.drawLine(cx, cy, yx, yy);
            g.drawLine(cx, cy, zx, zy);
            g.fillOval(xx-2, xy-2, 5, 5);
            g.fillOval(yx-2, yy-2, 5, 5);
            g.fillOval(zx-2, zy-2, 5, 5);
            g.drawString("x", xx-2, xy-4);
            g.drawString("y", yx-2, yy-4);
            g.drawString("z", zx-2, zy-4);
        }
    }


    public void advanceProps() {
        int[] pnum = this.animpropnum;

        for (int i = 0; i < pat.getNumberOfPaths(); i++)
            temppropnum[invpathperm.getMapping(i+1)-1] = pnum[i];
        for (int i = 0; i < pat.getNumberOfPaths(); i++)
            pnum[i] = temppropnum[i];
    }

    public void initAnimator() {
        findMaxMin();
        syncRenderersToSize();

        // figure out timing constants; adjust fps to get integer number of frames in loop
        num_frames = (int)(0.5 + (pat.getLoopEndTime() - pat.getLoopStartTime()) * jc.slowdown * jc.fps);
        sim_interval_secs = (pat.getLoopEndTime()-pat.getLoopStartTime()) / num_frames;
        real_interval_millis = (long)(1000.0 * sim_interval_secs * jc.slowdown);

        animpropnum = new int[pat.getNumberOfPaths()];
        for (int i = 1; i <= pat.getNumberOfPaths(); i++)
            animpropnum[i-1] = pat.getPropAssignment(i);
        temppropnum = new int[pat.getNumberOfPaths()];
        invpathperm = pat.getPathPermutation().getInverse();
    }

    private void findMaxMin() {
        // the algorithm here could be improved to take into account which props are
        // on which paths.  We may also want to leave room for the rest of the juggler.
        int i;
        Coordinate patternmax = null, patternmin = null;
        Coordinate handmax = null, handmin = null;
        Coordinate propmax = null, propmin = null;

        for (i = 1; i <= pat.getNumberOfPaths(); i++) {
            patternmax = Coordinate.max(patternmax, pat.getPathMax(i));
            patternmin = Coordinate.min(patternmin, pat.getPathMin(i));

            if (jugglinglab.core.Constants.DEBUG_LAYOUT)
                System.out.println("Pattern max "+i+" = "+patternmax);
        }

        // make sure all hands are visible
        for (i = 1; i <= pat.getNumberOfJugglers(); i++) {
            handmax = Coordinate.max(handmax, pat.getHandMax(i, HandLink.LEFT_HAND));
            handmin = Coordinate.min(handmin, pat.getHandMin(i, HandLink.LEFT_HAND));
            handmax = Coordinate.max(handmax, pat.getHandMax(i, HandLink.RIGHT_HAND));
            handmin = Coordinate.min(handmin, pat.getHandMin(i, HandLink.RIGHT_HAND));
        }

        for (i = 1; i <= pat.getNumberOfProps(); i++) {
            propmax = Coordinate.max(propmax, pat.getProp(i).getMax());
            propmin = Coordinate.min(propmin, pat.getProp(i).getMin());
        }

        // make sure props are entirely visible along all paths
        patternmax = Coordinate.add(patternmax, propmax);
        patternmin = Coordinate.add(patternmin, propmin);

        // make sure hands are entirely visible
        handmax = Coordinate.add(handmax, ren1.getHandWindowMax());
        handmin = Coordinate.add(handmin, ren1.getHandWindowMin());

        // make sure jugglers' bodies are visible
        this.overallmax = Coordinate.max(handmax, ren1.getJugglerWindowMax());
        this.overallmax = Coordinate.max(overallmax, patternmax);

        this.overallmin = Coordinate.min(handmin, ren1.getJugglerWindowMin());
        this.overallmin = Coordinate.min(overallmin, patternmin);

		// we want to ensure everything stays visible as we rotate the camera
		// viewpoint.  the following is simple and seems to work ok.
		if (pat.getNumberOfJugglers() == 1) {
			overallmin.z -= 0.3 * Math.max(Math.abs(overallmin.y), Math.abs(overallmax.y));
			overallmax.z += 5.0;	// keeps objects from rubbing against top of window
		} else {
			double tempx = Math.max(Math.abs(overallmin.x), Math.abs(overallmax.x));
			double tempy = Math.max(Math.abs(overallmin.y), Math.abs(overallmax.y));
			overallmin.z -= 0.4 * Math.max(tempx, tempy);
			overallmax.z += 0.4 * Math.max(tempx, tempy);
		}

		// make the x-coordinate origin at the center of the view
		double maxabsx = Math.max(Math.abs(this.overallmin.x), Math.abs(this.overallmax.x));
		this.overallmin.x = -maxabsx;
		this.overallmax.x = maxabsx;

        if (jugglinglab.core.Constants.DEBUG_LAYOUT) {
            System.out.println("Hand max = "+handmax);
            System.out.println("Hand min = "+handmin);
            System.out.println("Prop max = "+propmax);
            System.out.println("Prop min = "+propmin);
            System.out.println("Pattern max = "+patternmax);
            System.out.println("Pattern min = "+patternmin);
            System.out.println("Overall max = "+this.overallmax);
            System.out.println("Overall min = "+this.overallmin);

            this.overallmax = new Coordinate(100.0,0.0,500.0);
            this.overallmin = new Coordinate(-100.0,0.0,-100.0);
        }
    }

    private void syncRenderersToSize() {
        Dimension d = new Dimension(this.dim);

        if (this.jc.stereo) {
            d.width /= 2;
            this.ren1.initDisplay(d, jc.border, this.overallmax, this.overallmin);
            this.ren2.initDisplay(d, jc.border, this.overallmax, this.overallmin);
        } else
            this.ren1.initDisplay(d, jc.border, this.overallmax, this.overallmin);
    }


    public int[] getAnimPropNum() { return animpropnum; }

    public Color getBackground() { return ren1.getBackground(); }

    public AnimationPrefs getAnimationPrefs() { return jc; }


    // Helper method for writing animated GIFs
    // Adapted from https://community.oracle.com/thread/1264385
    public static void configureGIFMetadata(IIOMetadata meta,
                                            String delayTime,
                                            int imageIndex) {
        String metaFormat = meta.getNativeMetadataFormatName();

        if (!"javax_imageio_gif_image_1.0".equals(metaFormat)) {
            throw new IllegalArgumentException(
                    "Unfamiliar gif metadata format: " + metaFormat);
        }

        Node root = meta.getAsTree(metaFormat);

        //find the GraphicControlExtension node
        Node child = root.getFirstChild();
        while (child != null) {
            if ("GraphicControlExtension".equals(child.getNodeName())) {
                break;
            }
            child = child.getNextSibling();
        }

        IIOMetadataNode gce = (IIOMetadataNode) child;
        gce.setAttribute("userInputFlag", "FALSE");
        gce.setAttribute("delayTime", delayTime);

        //only the first node needs the ApplicationExtensions node
        if (imageIndex == 0) {
            IIOMetadataNode aes =
                    new IIOMetadataNode("ApplicationExtensions");
            IIOMetadataNode ae =
                    new IIOMetadataNode("ApplicationExtension");
            ae.setAttribute("applicationID", "NETSCAPE");
            ae.setAttribute("authenticationCode", "2.0");
            byte[] uo = new byte[]{
                //last two bytes is an unsigned short (little endian) that
                //indicates the the number of times to loop.
                //0 means loop forever.
                0x1, 0x0, 0x0
            };
            ae.setUserObject(uo);
            aes.appendChild(ae);
            root.appendChild(aes);
        }

        try {
            meta.setFromTree(metaFormat, root);
        } catch (IIOInvalidTreeException e) {
            //shouldn't happen
            throw new Error(e);
        }
    }


    // Called when the user wants to save a GIF from the command line. This skips
    // the file dialog box, progress monitor, and separate worker thread.
    public void writeGIF(OutputStream os, Animator.WriteGIFMonitor wgm) throws
                        IOException, JuggleExceptionInternal {

        ImageWriter iw = ImageIO.getImageWritersByFormatName("gif").next();
        ImageOutputStream ios = new MemoryCacheImageOutputStream(os);
        iw.setOutput(ios);
        iw.prepareWriteSequence(null);

        BufferedImage image = new BufferedImage(this.dim.width, this.dim.height,
                                                BufferedImage.TYPE_INT_RGB);
        Graphics g = image.getGraphics();

        int[] gifpropnum = new int[pat.getNumberOfPaths()];
        for (int i = 0; i < pat.getNumberOfPaths(); i++)
            gifpropnum[i] = pat.getPropAssignment(i+1);
        int patperiod = pat.getPeriod();
        int totalframes = patperiod * num_frames;
        int framecount = 0;

        String delayTime = String.valueOf((int)(real_interval_millis/10));

        for (int i = 0; i < patperiod; i++)  {
            double time = pat.getLoopStartTime();

            for (int j = 0; j < num_frames; j++) {
                this.drawFrame(time, g, false);

                ImageWriteParam iwp = iw.getDefaultWriteParam();
                IIOMetadata metadata = iw.getDefaultImageMetadata(
                        new ImageTypeSpecifier(image), iwp);
                configureGIFMetadata(metadata, delayTime, framecount);
                IIOImage ii = new IIOImage(image, null, metadata);
                iw.writeToSequence(ii, (ImageWriteParam) null);

                time += sim_interval_secs;
                framecount++;

                if (wgm != null) {
                    wgm.update(framecount, totalframes);
                    if (wgm.isCanceled()) {
                        ios.close();
                        os.close();
                        return;
                    }
                }
            }

            this.advanceProps();
        }

        g.dispose();
        iw.endWriteSequence();
        ios.close();
        os.close();
    }


    // Version that uses our own standalone GIF writer. It has trouble building
    // the color map when there are many individual colors, for example with the
    // image prop.
    /*
    public void writeGIF_old(OutputStream os, Animator.WriteGIFMonitor wgm) throws
                IOException, JuggleExceptionInternal {

        // Create the object that will actually do the writing
        GIFAnimWriter gaw = new GIFAnimWriter();

        int appWidth = this.dim.width;
        int appHeight = this.dim.height;

        BufferedImage image = new BufferedImage(appWidth, appHeight, BufferedImage.TYPE_INT_RGB);
        Graphics g = image.getGraphics();

        int[] gifpropnum = new int[pat.getNumberOfPaths()];
        for (int i = 0; i < pat.getNumberOfPaths(); i++)
            gifpropnum[i] = pat.getPropAssignment(i+1);
        int patperiod = pat.getPeriod();
        int totalframes = patperiod * num_frames * 2;
        int framecount = 0;

        // loop through the individual frames twice, first to build the
        // color map and the second to write the GIF frames
        for (int pass = 0; pass < 2; pass++) {
            if (pass == 1)
                gaw.writeHeader(os);

            for (int i = 0; i < patperiod; i++)  {
                double time = pat.getLoopStartTime();

                for (int j = 0; j < num_frames; j++) {
                    if (pass == 1)
                        gaw.writeDelay((int)(real_interval_millis/10), os);

                    this.drawFrame(time, g, false);

                    if (pass == 0)
                        gaw.doColorMap(image);
                    else
                        gaw.writeGIF(image, os);

                    if (wgm != null) {
                        framecount++;
                        wgm.update(framecount, totalframes);
                        if (wgm.isCanceled()) {
                            os.close();
                            return;
                        }
                    }

                    time += sim_interval_secs;
                }

                this.advanceProps();
            }
        }

        gaw.writeTrailer(os);
        g.dispose();
        os.close();
    }
    */

    public interface WriteGIFMonitor {
        // callback method invoked when a processing step is completed
        public void update(int step, int steps_total);

        // callback method returns true when user wants to cancel
        public boolean isCanceled();
    }
}
