// Juggler.java
//
// Copyright 2002-2022 Jack Boyce and the Juggling Lab contributors

package jugglinglab.renderer;

import jugglinglab.util.*;
import jugglinglab.jml.JMLPattern;
import jugglinglab.jml.HandLink;


// This class calculates the coordinates of the juggler elbows, shoulders, etc.

public class Juggler {
        // juggler dimensions, in centimeters
    public final static double shoulder_hw = 23.0;  // shoulder half-width (cm)
    public final static double shoulder_h = 40.0;  // throw pos. to shoulder
    public final static double waist_hw = 17.0;  // waist half-width
    public final static double waist_h = -5.0;
    public final static double elbow_hw = 30.0;  // elbow "home"
    public final static double elbow_h = 6.0;
    public final static double elbow_slop = 12.0;
    public final static double hand_out = 5.0;  // outside width of hand
    public final static double hand_in = 5.0;
    public final static double head_hw = 10.0;  // head half-width
    public final static double head_h = 26.0;  // head height
    public final static double neck_h = 5.0;  // neck height
    public final static double shoulder_y = 0.0;
    public final static double pattern_y = 30.0;
    public final static double upper_length = 41.0;
    public final static double lower_length = 40.0;

    public final static double lower_gap_wrist = 1.0;
    public final static double lower_gap_elbow = 0.0;
    public final static double lower_hand_height = 0.0;
    public final static double upper_gap_elbow = 0.0;
    public final static double upper_gap_shoulder = 0.0;

    protected final static double lower_total = lower_length + lower_gap_wrist + lower_gap_elbow;
    protected final static double upper_total = upper_length + upper_gap_elbow + upper_gap_shoulder;

        // the remaining are used only for the 3d display
    public final static double shoulder_radius = 6;
    public final static double elbow_radius = 4;
    public final static double wrist_radius = 2;


    public static void findJugglerCoordinates(JMLPattern pat, double time, JLVector[][] result)
                                    throws JuggleExceptionInternal {
        for (int juggler = 1; juggler <= pat.getNumberOfJugglers(); juggler++) {
            JLVector lefthand, righthand;
            JLVector leftshoulder, rightshoulder;
            JLVector leftelbow, rightelbow;
            JLVector leftwaist, rightwaist;
            JLVector leftheadbottom, leftheadtop;
            JLVector rightheadbottom, rightheadtop;

            Coordinate coord0 = new Coordinate();
            Coordinate coord1 = new Coordinate();
            Coordinate coord2 = new Coordinate();
            pat.getHandCoordinate(juggler, HandLink.LEFT_HAND, time, coord0);
            pat.getHandCoordinate(juggler, HandLink.RIGHT_HAND, time, coord1);
            lefthand = new JLVector(coord0.x,
                        coord0.z + lower_hand_height, coord0.y);
            righthand = new JLVector(coord1.x,
                        coord1.z + lower_hand_height, coord1.y);

            pat.getJugglerPosition(juggler, time, coord2);
            double angle = Math.toRadians(pat.getJugglerAngle(juggler, time));
            double s = Math.sin(angle);
            double c = Math.cos(angle);

            leftshoulder = new JLVector(
                coord2.x - shoulder_hw * c - shoulder_y * s,
                coord2.z + shoulder_h,
                coord2.y - shoulder_hw * s + shoulder_y * c);
            rightshoulder = new JLVector(
                coord2.x + shoulder_hw * c - shoulder_y * s,
                coord2.z + shoulder_h,
                coord2.y + shoulder_hw * s + shoulder_y * c);
            leftwaist = new JLVector(
                coord2.x - waist_hw * c - shoulder_y * s,
                coord2.z + waist_h,
                coord2.y - waist_hw * s + shoulder_y * c);
            rightwaist = new JLVector(
                coord2.x + waist_hw * c - shoulder_y * s,
                coord2.z + waist_h,
                coord2.y + waist_hw * s + shoulder_y * c);
            leftheadbottom = new JLVector(
                coord2.x - head_hw * c - shoulder_y * s,
                coord2.z + shoulder_h + neck_h,
                coord2.y - head_hw * s + shoulder_y * c);
            leftheadtop = new JLVector(
                coord2.x - head_hw * c - shoulder_y * s,
                coord2.z + shoulder_h + neck_h + head_h,
                coord2.y - head_hw * s + shoulder_y * c);
            rightheadbottom = new JLVector(
                coord2.x + head_hw * c - shoulder_y * s,
                coord2.z + shoulder_h + neck_h,
                coord2.y + head_hw * s + shoulder_y * c);
            rightheadtop = new JLVector(
                coord2.x + head_hw * c - shoulder_y * s,
                coord2.z + shoulder_h + neck_h + head_h,
                coord2.y + head_hw * s + shoulder_y * c);

            double L = lower_total; // length of the lower arm
            double U = upper_total; // length of the upper arm
            JLVector deltaL = JLVector.sub(lefthand, leftshoulder);
            double D = deltaL.length();
            if (D <= (L+U)) {
                // Calculate the coordinates of the elbows
                double Lr = Math.sqrt((4.0*U*U*L*L-(U*U+L*L-D*D)*(U*U+L*L-D*D))/(4.0*D*D));
                if (Double.isNaN(Lr))
                    throw new JuggleExceptionInternal("NaN in renderer 1");

                double factor = Math.sqrt(U*U-Lr*Lr)/D;
                if (Double.isNaN(factor))
                    throw new JuggleExceptionInternal("NaN in renderer 2");
                JLVector Lxsc = JLVector.scale(factor, deltaL);
                double Lalpha = Math.asin(deltaL.y / D);
                if (Double.isNaN(Lalpha))
                    throw new JuggleExceptionInternal("NaN in renderer 3");
                factor = 1.0 + Lr*Math.tan(Lalpha)/(factor*D);
                leftelbow = new JLVector(
                        leftshoulder.x + Lxsc.x * factor,
                        leftshoulder.y + Lxsc.y - Lr*Math.cos(Lalpha),
                        leftshoulder.z + Lxsc.z * factor);
            } else {
                // if the distance between hand and shoulder is unrealistically far, the elbow is halfway between them
                leftelbow = new JLVector(
                    (leftshoulder.x + lefthand.x) / 2,
                    (leftshoulder.y + lefthand.y) / 2,
                    (leftshoulder.z + lefthand.z) / 2);
            }

            JLVector deltaR = JLVector.sub(righthand, rightshoulder);
            D = deltaR.length();
            if (D <= (L+U)) {
                // Calculate the coordinates of the elbows
                double Rr = Math.sqrt((4.0*U*U*L*L-(U*U+L*L-D*D)*(U*U+L*L-D*D))/(4.0*D*D));
                if (Double.isNaN(Rr))
                    throw new JuggleExceptionInternal("NaN in renderer 4");

                double factor = Math.sqrt(U*U-Rr*Rr)/D;
                if (Double.isNaN(factor))
                    throw new JuggleExceptionInternal("NaN in renderer 5");
                JLVector Rxsc = JLVector.scale(factor, deltaR);
                double Ralpha = Math.asin(deltaR.y / D);
                if (Double.isNaN(Ralpha))
                    throw new JuggleExceptionInternal("NaN in renderer 6");
                factor = 1.0 + Rr*Math.tan(Ralpha)/(factor*D);
                rightelbow = new JLVector(
                        rightshoulder.x + Rxsc.x * factor,
                        rightshoulder.y + Rxsc.y - Rr*Math.cos(Ralpha),
                        rightshoulder.z + Rxsc.z * factor);
            } else {
                // if the distance between hand and shoulder is unrealistically far, the elbow is halfway between them
                rightelbow = new JLVector(
                    (rightshoulder.x + righthand.x) / 2,
                    (rightshoulder.y + righthand.y) / 2,
                    (rightshoulder.z + righthand.z) / 2);
            }

            result[juggler - 1][0] = lefthand;
            result[juggler - 1][1] = righthand;
            result[juggler - 1][2] = leftshoulder;
            result[juggler - 1][3] = rightshoulder;
            result[juggler - 1][4] = leftelbow;
            result[juggler - 1][5] = rightelbow;
            result[juggler - 1][6] = leftwaist;
            result[juggler - 1][7] = rightwaist;
            result[juggler - 1][8] = leftheadbottom;
            result[juggler - 1][9] = leftheadtop;
            result[juggler - 1][10] = rightheadbottom;
            result[juggler - 1][11] = rightheadtop;
        }
    }

}
