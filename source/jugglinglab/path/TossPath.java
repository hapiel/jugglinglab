// TossPath.java
//
// Copyright 2019 by Jack Boyce (jboyce@gmail.com)

package jugglinglab.path;

import java.text.MessageFormat;
import jugglinglab.util.*;


public class TossPath extends Path {
    protected static final double g_def = 980;  // using CGS units

    protected double    bx, cx;
    protected double    by, cy;
    protected double    az, bz, cz;

    protected double    g = g_def;

    @Override
    public String getType() { return "Toss"; }

    @Override
    public ParameterDescriptor[] getParameterDescriptors() {
        ParameterDescriptor[] result = new ParameterDescriptor[1];

        result[0] = new ParameterDescriptor("g", ParameterDescriptor.TYPE_FLOAT,
                                            null, new Double(g_def), new Double(g));
        return result;
    }

    @Override
    public void initPath(String st) throws JuggleExceptionUser {
        double g = g_def;

        // now parse for edits to the above variables
        ParameterList pl = new ParameterList(st);
        for (int i = 0; i < pl.getNumberOfParameters(); i++) {
            String pname = pl.getParameterName(i);
            String pvalue = pl.getParameterValue(i);

            if (pname.equalsIgnoreCase("g")) {
                try {
                    g = Double.valueOf(pvalue).doubleValue();
                } catch (NumberFormatException nfe) {
                    String template = errorstrings.getString("Error_number_format");
                    Object[] arguments = { "g" };
                    throw new JuggleExceptionUser(MessageFormat.format(template, arguments));
                }
            } else
                throw new JuggleExceptionUser(errorstrings.getString("Error_path_badmod")+": '"+pname+"'");
        }
        this.g = g;
        az = -0.5 * g;
    }

    @Override
    public void calcPath() throws JuggleExceptionInternal {
        if (start_coord == null || end_coord == null)
            throw new JuggleExceptionInternal("Error in parabolic path: endpoints not set");

        double  t = getDuration();
        cx = start_coord.x;
        bx = (end_coord.x - start_coord.x) / t;
        cy = start_coord.y;
        by = (end_coord.y - start_coord.y) / t;
        cz = start_coord.z;
        bz = (end_coord.z - start_coord.z) / t - az * t;
    }

    @Override
    public Coordinate getStartVelocity() {
        return new Coordinate(bx, by, bz);
    }

    @Override
    public Coordinate getEndVelocity() {
        return new Coordinate(bx, by, bz+2.0*az*getDuration());
    }

    @Override
    public void getCoordinate(double time, Coordinate newPosition) {
        if ((time < start_time) || (time > end_time))
            return;
        time -= start_time;
        newPosition.setCoordinate(cx+bx*time, cy+by*time, cz+time*(bz+az*time));
    }

    @Override
    protected Coordinate getMax2(double begin, double end) {
        Coordinate result = null;
        double tlow = Math.max(start_time, begin);
        double thigh = Math.min(end_time, end);

        result = check(result, tlow, true);
        result = check(result, thigh, true);

        if (az < 0.0) {
            double te = -bz / (2.0*az) + start_time;
            if (tlow < te && te < thigh)
                result = check(result, te, true);
        }
        return result;
    }

    @Override
    protected Coordinate getMin2(double begin, double end) {
        Coordinate result = null;
        double tlow = Math.max(start_time, begin);
        double thigh = Math.min(end_time, end);

        result = check(result, tlow, false);
        result = check(result, thigh, false);

        if (az > 0.0) {
            double te = -by / (2.0*az) + start_time;
            if (tlow < te && te < thigh)
                result = check(result, te, false);
        }
        return result;
    }
}