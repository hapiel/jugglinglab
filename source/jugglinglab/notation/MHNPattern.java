// MHNPattern.java
//
// Copyright 2020 by Jack Boyce (jboyce@gmail.com)

package jugglinglab.notation;

import java.text.MessageFormat;
import java.util.*;
import java.awt.event.*;

import jugglinglab.core.*;
import jugglinglab.util.*;
import jugglinglab.jml.*;


// "Multi-Hand Notation" (name due to Ed Carstens) is a very general notation that
// defines a matrix of hands and discrete event times, and describes patterns as
// transitions between elements in this matrix.
//
// MHN has no standardized string (textual) representation. Hence this class is
// abstract because it lacks a fromString() method. It is however useful as a
// building block for other notations. For example we model siteswap notation
// as a type of MHN, with a parser to interpret siteswap notation into the
// internal MHN matrix representation.
//
// The main business of this class is to go from the matrix-based internal
// representation of an MHN pattern to JML. Any notation that builds on MHN can
// avoid duplicating this functionality.

public abstract class MHNPattern extends Pattern {
    protected static double bps_default = -1.0;  // calculate bps
    protected static double dwell_default = 1.3;
    protected static double gravity_default = 980.0;
    protected static double propdiam_default = 10.0;
    protected static double bouncefrac_default = 0.9;
    protected static String prop_default = "ball";
    //Mahit begin
    protected static boolean hold_default = false;
    protected static boolean dwellmax_default = true;
    //Mahit end

    // original config string:
    protected String config;

    // input parameters:
    protected String pattern;
    protected double bps = bps_default;
    protected double dwell = dwell_default;
    protected double gravity = gravity_default;
    protected double propdiam = propdiam_default;
    protected double bouncefrac = bouncefrac_default;
    protected String prop = prop_default;
    protected String[] color;
    protected String title;
    //Mahit begin
    protected String hss;  
    protected boolean hold = hold_default; 
    protected boolean dwellmax = dwellmax_default;
    protected String handspec;
    protected double[] dwellarray;
    //Mahit end

    // internal variables:
    protected int numjugglers;
    protected int numpaths;
    protected int period;
    protected int max_occupancy;
    protected MHNThrow[][][][] th;
    protected MHNHands hands;
    protected MHNBody bodies;
    protected int max_throw;
    protected int indexes;
    protected ArrayList<MHNSymmetry> symmetry;

    public static final int RIGHT_HAND = 0;
    public static final int LEFT_HAND = 1;

    public int getNumberOfJugglers()         { return numjugglers; }
    public int getNumberOfPaths()            { return numpaths; }
    public int getPeriod()                   { return period; }
    public int getIndexes()                  { return indexes; }
    public int getMaxOccupancy()             { return max_occupancy; }
    public int getMaxThrow()                 { return max_throw; }
    public MHNThrow[][][][] getThrows()      { return th; }
    public int getNumberOfSymmetries()       { return symmetry.size(); }
    public String getPropName()              { return prop; }
    public void addSymmetry(MHNSymmetry ss)  { symmetry.add(ss); }
    public MHNSymmetry getSymmetry(int i)    { return symmetry.get(i); }


    // pull out the MHN-related parameters from the given list, leaving any
    // other parameters alone.
    @Override
    public MHNPattern fromParameters(ParameterList pl) throws
                                    JuggleExceptionUser, JuggleExceptionInternal {
        this.config = pl.toString();                // save for toString()

        pattern = pl.removeParameter("pattern");    // only required parameter
        if (pattern == null)
            throw new JuggleExceptionUser(errorstrings.getString("Error_no_pattern"));

        String temp = null;
        if ((temp = pl.removeParameter("bps")) != null) {
            try {
                bps = Double.parseDouble(temp);
            } catch (NumberFormatException nfe) {
                throw new JuggleExceptionUser(errorstrings.getString("Error_bps_value"));
            }
        }

        if ((temp = pl.removeParameter("dwell")) != null) {
            try {
                dwell = Double.parseDouble(temp);
            } catch (NumberFormatException nfe) {
                throw new JuggleExceptionUser(errorstrings.getString("Error_dwell_value"));
            }
        }

        if ((temp = pl.removeParameter("hands")) != null)
            hands = new MHNHands(temp);

        if ((temp = pl.removeParameter("body")) != null)
            bodies = new MHNBody(temp);

        if ((temp = pl.removeParameter("gravity")) != null) {
            try {
                gravity = Double.parseDouble(temp);
            } catch (NumberFormatException e) {
            }
        }

        if ((temp = pl.removeParameter("propdiam")) != null) {
            try {
                propdiam = Double.parseDouble(temp);
            } catch (NumberFormatException e) {
            }
        }

        if ((temp = pl.removeParameter("bouncefrac")) != null) {
            try {
                bouncefrac = Double.parseDouble(temp);
            } catch (NumberFormatException e) {
            }
        }

        if ((temp = pl.removeParameter("prop")) != null) {
            prop = temp;
        }

        if ((temp = pl.removeParameter("colors")) != null) {
            if (temp.trim().equals("mixed"))
                temp = "{red}{green}{blue}{yellow}{cyan}{magenta}{orange}{pink}{gray}{black}";
            else
                temp = JLFunc.expandRepeats(temp);

            StringTokenizer st1 = new StringTokenizer(temp, "}", false);
            StringTokenizer st2 = null;
            String          str = null;

            int numcolors = st1.countTokens();
            color = new String[numcolors];

            // Parse the colors parameter
            for (int i = 0; i < numcolors; i++) {
                // Look for next {...} block
                str = st1.nextToken().replace('{', ' ').trim();

                // Parse commas
                st2 = new StringTokenizer(str, ",", false);

                switch (st2.countTokens()) {
                    case 1:
                        // Use the value as a color name
                        color[i] = st2.nextToken().trim().toLowerCase();
                        break;
                    case 3:
                        // Use the three values as RGB values
                        color[i] = "{" + str + "}";
                        break;
                    default:
                        throw new JuggleExceptionUser(errorstrings.getString("Error_color_format"));
                }
            }
        }

//Mahit begin
        if ((temp = pl.removeParameter("hss")) != null) {
            hss = temp;
        }

        if (hss != null) {
            if ((temp = pl.removeParameter("hold")) != null) {
                try {
                    hold = Boolean.parseBoolean(temp);
                } catch (IllegalFormatException ife) {
                    throw new JuggleExceptionUser("Please enter true/false for hold");
                }
            } 

            if ((temp = pl.removeParameter("dwellmax")) != null) {
                try {
                    dwellmax = Boolean.parseBoolean(temp);
                } catch (IllegalFormatException ife) {
                    throw new JuggleExceptionUser("Please enter true/false for dwellmax");
                }
            }

            if ((temp = pl.removeParameter("handspec")) != null) {
                handspec = temp;
            }
        }
//Mahit end
        if ((temp = pl.removeParameter("title")) != null) {
            this.title = temp.trim();
        }

        return this;
    }

    @Override
    public String toString() {
        // print out configuration parameters in a standard order
        if (this.config == null)
            return null;

        ParameterList pl = new ParameterList(this.config);
        String result = "";

        // write the parameters out in a standard order
        List<String> keys = Arrays.asList("pattern", "bps", "dwell", "hands", "body",
                                "gravity", "propdiam", "bouncefrac", "prop", "colors",
                                "title");

        for (String key : keys) {
            String value = pl.getParameter(key);
            if (value != null)
                result += key + "=" + value + ";";
        }

        if (result.length() > 0)
            result = result.substring(0, result.length() - 1);

        return result;
    }

    //--------------------------------------------------------------------------
    // Build out internal representation from parsed pattern
    //--------------------------------------------------------------------------

    protected void buildRepresentation() throws JuggleExceptionUser, JuggleExceptionInternal {
        // build out the internal pattern representation in steps
        //
        // this will find and raise any errors in the pattern
        if (Constants.DEBUG_PARSING) {
            System.out.println("-----------------------------------------------------");
            System.out.println("Building internal MHNPattern representation...\n");
            System.out.println("findMasterThrows()");
        }
        findMasterThrows();

        if (Constants.DEBUG_PARSING)
            System.out.println("assignPaths()");
        assignPaths();

        if (Constants.DEBUG_PARSING)
            System.out.println("findThrowSources()");
        findThrowSources();

        if (Constants.DEBUG_PARSING)
            System.out.println("setCatchOrder()");
        setCatchOrder();

        if (Constants.DEBUG_PARSING) {
            String s = getInternalRepresentation();
            if (s != null) {
                System.out.println("\nInternal MHNPattern representation:\n");
                System.out.println(s);
                System.out.println("-----------------------------------------------------");
            }
        }
    }

    // Determines which throws are "master" throws. Because of symmetries defined
    // for the pattern, some throws are shifted or reflected copies of others.
    // For each such chain of related throws, appoint one as the master.
    protected void findMasterThrows() throws JuggleExceptionInternal {
        // start by making every throw a master throw
        for (int i = 0; i < indexes; ++i) {
            for (int j = 0; j < numjugglers; ++j) {
                for (int h = 0; h < 2; ++h) {
                    for (int slot = 0; slot < max_occupancy; ++slot) {
                        MHNThrow mhnt = th[j][h][i][slot];
                        if (mhnt != null) {
                            mhnt.master = mhnt;
                            mhnt.source = null;
                        }
                    }
                }
            }
        }

        boolean changed = true;
        while (changed) {
            changed = false;

            for (int s = 0; s < getNumberOfSymmetries(); ++s) {
                MHNSymmetry sym = getSymmetry(s);
                Permutation jperm = sym.getJugglerPerm();
                int delay = sym.getDelay();

                for (int i = 0; i < indexes; ++i) {
                    int imagei = i + delay;
                    if (imagei >= indexes)
                        continue;

                    for (int j = 0; j < numjugglers; ++j) {
                        for (int h = 0; h < 2; ++h) {
                            for (int slot = 0; slot < max_occupancy; ++slot) {
                                MHNThrow mhnt = th[j][h][i][slot];
                                if (mhnt == null)
                                    continue;

                                int imagej = jperm.getMapping(j + 1);
                                int imageh = (imagej > 0 ? h : 1 - h);
                                imagej = Math.abs(imagej) - 1;

                                MHNThrow imaget = th[imagej][imageh][imagei][slot];
                                if (imaget == null)
                                    throw new JuggleExceptionInternal("Problem finding master throws");

                                MHNThrow m = mhnt.master;
                                MHNThrow im = imaget.master;
                                if (m == im)
                                    continue;

                                // we have a disagreement about which is the master;
                                // choose one of them and set them equal
                                MHNThrow newm = m;
                                if (m.index > im.index)
                                    newm = im;
                                else if (m.index == im.index) {
                                    if (m.juggler > im.juggler)
                                        newm = im;
                                    else if (m.juggler == im.juggler) {
                                        if (m.hand > im.hand)
                                            newm = im;
                                    }
                                }
                                mhnt.master = imaget.master = newm;
                                changed = true;
                            }
                        }
                    }
                }
            }
        }

        if (Constants.DEBUG_PARSING) {
            for (int i = 0; i < indexes; ++i) {
                for (int j = 0; j < numjugglers; ++j) {
                    for (int h = 0; h < 2; ++h) {
                        for (int slot = 0; slot < max_occupancy; ++slot) {
                            MHNThrow mhnt = th[j][h][i][slot];
                            if (mhnt != null && mhnt.master == mhnt)
                                System.out.println("master throw at j="+j+",h="+h+",i="+i+",slot="+slot);
                        }
                    }
                }
            }
        }
    }

    // Figures out which throws are filling other throws, and assigns path
    // numbers to all throws.
    //
    // This process is complicated by the fact that we allow multiplexing.
    protected void assignPaths() throws JuggleExceptionUser, JuggleExceptionInternal {
        for (int i = 0; i < indexes; ++i) {
            for (int j = 0; j < numjugglers; ++j) {
                for (int h = 0; h < 2; ++h) {
                    for (int slot = 0; slot < max_occupancy; ++slot) {
                        MHNThrow sst = th[j][h][i][slot];
                        if (sst == null || sst.master != sst)  // loop over master throws
                            continue;

                        // Figure out which slot number we're filling with this
                        // master throw. We need to find a value of `targetslot`
                        // that is empty for the master throw's target, as well as
                        // the targets of its images.

                        int targetslot = 0;
                        while (targetslot < max_occupancy) {  // find value of targetslot that works
                            boolean itworks = true;

                            // loop over all throws that have sst as master
                            for (int i2 = 0; i2 < indexes; ++i2) {
                                for (int j2 = 0; j2 < numjugglers; ++j2) {
                                    for (int h2 = 0; h2 < 2; ++h2) {
                                        for (int slot2 = 0; slot2 < max_occupancy; ++slot2) {
                                            MHNThrow sst2 = th[j2][h2][i2][slot2];
                                            if (sst2 == null || sst2.master != sst || sst2.targetindex >= indexes)
                                                continue;

                                            MHNThrow target = th[sst2.targetjuggler-1][sst2.targethand][sst2.targetindex][targetslot];
                                            if (target == null)
                                                itworks = false;
                                            else
                                                itworks &= (target.source == null);  // target also unfilled?
                                        }
                                    }
                                }
                            }
                            if (itworks)
                                break;
                            ++targetslot;
                        }

                        if (targetslot == max_occupancy) {
                            if (Constants.DEBUG_PARSING)
                                System.out.println("Error: Too many objects landing on beat "
                                        + (sst.targetindex + 1) + " for juggler "
                                        + sst.targetjuggler + ", "
                                        + (sst.targethand == 0 ? "right hand" : "left hand"));

                            String template = errorstrings.getString("Error_badpattern_landings");
                            String hand = (sst.targethand == 0 ? errorstrings.getString("Error_right_hand")
                                        : errorstrings.getString("Error_left_hand"));
                            Object[] arguments = { new Integer(sst.targetindex + 1),
                                        new Integer(sst.targetjuggler), hand };
                            throw new JuggleExceptionUser(MessageFormat.format(template, arguments));
                        }

                        // loop again over all throws that have sst as master,
                        // wiring up sources and targets using the value of `targetslot`

                        for (int i2 = 0; i2 < indexes; ++i2) {
                            for (int j2 = 0; j2 < numjugglers; ++j2) {
                                for (int h2 = 0; h2 < 2; ++h2) {
                                    for (int slot2 = 0; slot2 < max_occupancy; ++slot2) {
                                        MHNThrow sst2 = th[j2][h2][i2][slot2];
                                        if (sst2 == null || sst2.master != sst || sst2.targetindex >= indexes)
                                            continue;

                                        MHNThrow target2 = th[sst2.targetjuggler-1][sst2.targethand][sst2.targetindex][targetslot];
                                        if (target2 == null)
                                            throw new JuggleExceptionInternal("Got null target in assignPaths()");

                                        sst2.target = target2;      // hook source and target together
                                        target2.source = sst2;
                                    }
                                }
                            }
                        }


                    }
                }
            }
        }

        // assign path numbers to all of the throws
        int currentpath = 1;
        for (int i = 0; i < indexes; ++i) {
            for (int j = 0; j < numjugglers; ++j) {
                for (int h = 0; h < 2; ++h) {
                    for (int slot = 0; slot < max_occupancy; ++slot) {
                        MHNThrow sst = th[j][h][i][slot];
                        if (sst == null)
                            continue;

                        if (sst.source != null) {
                            sst.pathnum = sst.source.pathnum;
                            continue;
                        }

                        if (currentpath > numpaths) {
                            if (Constants.DEBUG_PARSING) {
                                System.out.println("j="+j+", h="+h+", index="+i+", slot="+slot+"\n");
                                System.out.println("---------------------------");
                                for (int tempi = 0; tempi <= i; ++tempi) {
                                    for (int temph = 0; temph < 2; ++temph) {
                                        for (int tempslot = 0; tempslot < max_occupancy; ++tempslot) {
                                            MHNThrow tempsst = th[0][temph][tempi][tempslot];
                                            System.out.println("index="+tempi+", hand="+temph+", slot="+tempslot+":");
                                            if (tempsst == null)
                                                System.out.println("   null entry");
                                            else
                                                System.out.println("   targetindex="+tempsst.targetindex+", targethand="+
                                                                   tempsst.targethand+", targetslot="+tempsst.targetslot+
                                                                   ", pathnum="+tempsst.pathnum);
                                        }
                                    }
                                }
                                System.out.println("---------------------------");
                            }

                            throw new JuggleExceptionUser(errorstrings.getString("Error_badpattern"));
                        }

                        sst.pathnum = currentpath;
                        ++currentpath;
                    }
                }
            }
        }

        if (currentpath <= numpaths)
            throw new JuggleExceptionInternal("Problem assigning path numbers 2");
    }

    // Sets the `source` field for all throws that don't already have it set.
    protected void findThrowSources() throws JuggleExceptionInternal {
        for (int i = indexes - 1; i >= 0; --i) {
            for (int j = 0; j < numjugglers; ++j) {
                for (int h = 0; h < 2; ++h) {
                    for (int slot = 0; slot < max_occupancy; ++slot) {
                        MHNThrow sst = th[j][h][i][slot];
                        if (sst == null || sst.source != null)
                            continue;

                        if ((i + getPeriod()) >= indexes)
                            throw new JuggleExceptionInternal("Could not get throw source 2");

                        MHNThrow sst2 = th[j][h][i + getPeriod()][slot].source;
                        if (sst2 == null)
                            throw new JuggleExceptionInternal("Could not get throw source 1");

                        MHNThrow sst3 = new MHNThrow();
                        sst3.juggler = sst2.juggler;
                        sst3.hand = sst2.hand;
                        sst3.index = sst2.index - getPeriod();
                        sst3.slot = sst2.slot;
                        sst3.targetjuggler = j;
                        sst3.targethand = h;
                        sst3.targetindex = i;
                        sst3.targetslot = slot;
                        sst3.handsindex = -1;   // undefined
                        sst3.pathnum = sst.pathnum;
                        sst3.mod = sst2.mod;
                        sst3.master = sst2.master;
                        sst3.source = null;
                        sst3.target = sst;

                        sst.source = sst3;
                    }
                }
            }
        }
    }

    // Decide whether the catches immediately prior to the two given throws should be
    // made in the order given, or whether they should be switched.
    //
    // The following implementation isn't ideal; we would like a function that is
    // invariant with respect to the various pattern symmetries we can apply, but
    // I don't think this is possible with respect to the jugglers.
    protected static boolean isCatchOrderIncorrect(MHNThrow t1, MHNThrow t2) {
        // first look at the time spent in the air; catch higher throws first
        if (t1.source.index > t2.source.index)
            return true;
        if (t1.source.index < t2.source.index)
            return false;

        // look at which juggler it's from; catch from "faraway" jugglers first
        int jdiff1 = Math.abs(t1.juggler - t1.source.juggler);
        int jdiff2 = Math.abs(t2.juggler - t2.source.juggler);
        if (jdiff1 < jdiff2)
            return true;
        if (jdiff1 > jdiff2)
            return false;

        // look at which hand it's from; catch from same hand first
        int hdiff1 = Math.abs(t1.hand - t1.source.hand);
        int hdiff2 = Math.abs(t2.hand - t2.source.hand);
        if (hdiff1 > hdiff2)
            return true;
        if (hdiff1 < hdiff2)
            return false;

        return false;
    }

    // Sets the MHNThrow.catching and MHNThrow.catchnum fields
    protected void setCatchOrder() throws JuggleExceptionInternal {
        MHNThrow[][][][] th = getThrows();

        // figure out the correct catch order for master throws...
        for (int k = 0; k < getIndexes(); ++k) {
            for (int j = 0; j < getNumberOfJugglers(); ++j) {
                for (int h = 0; h < 2; ++h) {
                    int slotcatches = 0;

                    for (int slot = 0; slot < getMaxOccupancy(); ++slot) {
                        MHNThrow sst = th[j][h][k][slot];
                        if (sst == null)
                            break;

                        sst.catching = (sst.source.mod.charAt(0) != 'H');
                        if (sst.catching) {
                            sst.catchnum = slotcatches;
                            ++slotcatches;
                        }
                    }

                    // Arrange the order of the catches, if more than one

                    if (slotcatches < 2)
                        continue;

                    for (int slot1 = 0; slot1 < getMaxOccupancy(); ++slot1) {
                        MHNThrow sst1 = th[j][h][k][slot1];
                        if (sst1 == null || sst1.master != sst1)  // only master throws
                            break;

                        if (!sst1.catching)
                            continue;

                        for (int slot2 = (slot1 + 1); slot2 < getMaxOccupancy(); ++slot2) {
                            MHNThrow sst2 = th[j][h][k][slot2];
                            if (sst2 == null || sst2.master != sst2)
                                break;

                            if (!sst2.catching)
                                continue;

                            boolean switchcatches = false;
                            if (sst1.catchnum < sst2.catchnum)
                                switchcatches = isCatchOrderIncorrect(sst1, sst2);
                            else
                                switchcatches = isCatchOrderIncorrect(sst2, sst1);

                            if (switchcatches) {
                                int temp = sst1.catchnum;
                                sst1.catchnum = sst2.catchnum;
                                sst2.catchnum = temp;
                            }
                        }
                    }
                }
            }
        }

        // ...and then copy that over to the non-master throws
        for (int k = 0; k < getIndexes(); ++k) {
            for (int j = 0; j < getNumberOfJugglers(); ++j) {
                for (int h = 0; h < 2; ++h) {
                    for (int slot = 0; slot < getMaxOccupancy(); ++slot) {
                        MHNThrow sst = th[j][h][k][slot];
                        if (sst == null || sst.master == sst)  // skip master throws
                            break;

                        sst.catchnum = sst.master.catchnum;
                    }
                }
            }
        }
    }

    // dump the internal state of the pattern to a string; intended to be used
    // for debugging
    protected String getInternalRepresentation() {
        StringBuffer sb = new StringBuffer();

        sb.append("numjugglers = " + getNumberOfJugglers() + "\n");
        sb.append("numpaths = " + getNumberOfPaths() + "\n");
        sb.append("period = " + getPeriod() + "\n");
        sb.append("max_occupancy = " + getMaxOccupancy() + "\n");
        sb.append("max_throw = " + getMaxThrow() + "\n");
        sb.append("indexes = " + getIndexes() + "\n");
        sb.append("throws:\n");

        for (int i = 0; i < getIndexes(); ++i) {
            for (int j = 0; j < numjugglers; ++j) {
                for (int h = 0; h < 2; ++h) {
                    for (int s = 0; s < getMaxOccupancy(); ++s) {
                        sb.append("  th[" + j + "][" + h + "][" + i + "][" + s + "] = ");
                        MHNThrow mhnt = th[j][h][i][s];
                        if (mhnt == null)
                            sb.append("null\n");
                        else
                            sb.append(mhnt.toString() + "\n");
                    }
                }
            }
        }

        sb.append("symmetries:\n");  // not finished
        sb.append("hands:\n");
        sb.append("bodies:\n");

        return sb.toString();
    }

    // returns the pattern's starting state
    //
    // result is a matrix of dimension (jugglers) x 2 x (indexes), with values
    // between 0 and (max_occupancy) inclusive
    public int[][][] getStartingState(int indexes) {
        int[][][] result = new int[getNumberOfJugglers()][2][indexes];

        for (int i = getPeriod(); i < getIndexes(); ++i) {
            for (int j = 0; j < numjugglers; ++j) {
                for (int h = 0; h < 2; ++h) {
                    for (int s = 0; s < getMaxOccupancy(); ++s) {
                        MHNThrow mhnt = th[j][h][i][s];

                        if (mhnt != null && mhnt.source.index < getPeriod()) {
                            if ((i - getPeriod()) < indexes)
                                result[j][h][i - getPeriod()]++;
                        }
                    }
                }
            }
        }

        return result;
    }

    //--------------------------------------------------------------------------
    // Convert from internal representation to JML
    //--------------------------------------------------------------------------

    // The following are default spatial coordinates to use
    protected static final double[] samethrowx =
        { 0.0, 20.0, 25.0, 12.0, 10.0,  7.5,  5.0,  5.0,  5.0 };
    protected static final double[] crossingthrowx =
        { 0.0, 20.0, 25.0, 12.0, 10.0, 18.0, 25.0, 25.0, 30.0 };
    protected static final double[] catchx =
        { 0.0, 30.0, 25.0, 30.0, 40.0, 45.0, 45.0, 50.0, 50.0 };
    protected static final double restingx = 25.0;

    // The following is used when multiple catches are made in a hand, on the same beat.
    // It specifies the fraction of a beat to spread the catches over.
    protected static final double splitcatchfactor = 0.4;

    @Override
    public JMLPattern asJMLPattern() throws JuggleExceptionUser, JuggleExceptionInternal {
        JMLPattern result = new JMLPattern();

        // Step 1 -- Set up the basic pattern information:

        // pattern title needs to be set elsewhere

        result.setNumberOfJugglers(getNumberOfJugglers());
        result.setNumberOfPaths(getNumberOfPaths());

        if (bps <= 0.0)       // signal that we should calculate bps
            bps = calcBps();

//Mahit begin
//        if (hss != null) {
//        	bps *= result.getNumberOfJugglers();
//        }
//Mahit end
        
        int balls = getNumberOfPaths();
        int props = (color == null) ? 1 : Math.min(balls, color.length);
        for (int i = 0; i < props; i++) {
            String mod = null;
            if (propdiam != MHNPattern.propdiam_default)
                mod = "diam="+propdiam;
            if (color != null) {
                String colorstr = "color="+color[i];
                if (mod == null)
                    mod = colorstr;
                else
                    mod = mod + ";" + colorstr;
            }
            result.addProp(new PropDef(getPropName(), mod));
        }
        int[] pa = new int[balls];
        for (int i = 0; i < balls; i++)
            pa[i] = 1 + (i % props);
        result.setPropAssignments(pa);

        // Step 2 -- Add the symmetries to the pattern:

        for (int i = 0; i < getNumberOfSymmetries(); i++) {
            MHNSymmetry sss = getSymmetry(i);
            int symtype;
            int[] pathmap = new int[balls+1];

            switch (sss.getType()) {
                case MHNSymmetry.TYPE_DELAY:
                    symtype = JMLSymmetry.TYPE_DELAY;
                    // figure out the path mapping
                    {
                        MHNThrow[][][][] th = getThrows();
                        for (int k = 0; k < (getIndexes()-sss.getDelay()); k++) {
                            for (int j = 0; j < getNumberOfJugglers(); j++) {
                                for (int h = 0; h < 2; h++) {
                                    for (int slot = 0; slot < getMaxOccupancy(); slot++) {
                                        MHNThrow sst = th[j][h][k][slot];
                                        if (sst != null && sst.pathnum != -1) {
                                            MHNThrow sst2 = th[j][h][k+sss.getDelay()][slot];
                                            if (sst2 == null)
                                                throw new JuggleExceptionUser(errorstrings.getString("Error_badpattern_paths"));
                                            if ((sst.pathnum == 0) || (sst2.pathnum == 0))
                                                throw new JuggleExceptionUser(errorstrings.getString("Error_badpattern_paths"));
                                            if (pathmap[sst.pathnum] == 0)
                                                pathmap[sst.pathnum] = sst2.pathnum;
                                            else if (pathmap[sst.pathnum] != sst2.pathnum)
                                                throw new JuggleExceptionUser(errorstrings.getString("Error_badpattern_delay"));
                                        }
                                    }
                                }
                            }
                        }
                    }
                        break;
                case MHNSymmetry.TYPE_SWITCH:
                    symtype = JMLSymmetry.TYPE_SWITCH;
                    break;
                case MHNSymmetry.TYPE_SWITCHDELAY:
                    symtype = JMLSymmetry.TYPE_SWITCHDELAY;

                    // figure out the path mapping
                    {
                        Permutation jugperm = sss.getJugglerPerm();

                        MHNThrow[][][][] th = getThrows();
                        for (int k = 0; k < (getIndexes()-sss.getDelay()); k++) {
                            for (int j = 0; j < getNumberOfJugglers(); j++) {
                                for (int h = 0; h < 2; h++) {
                                    for (int slot = 0; slot < getMaxOccupancy(); slot++) {
                                        MHNThrow sst = th[j][h][k][slot];
                                        if (sst != null && sst.pathnum != -1) {
                                            int map = jugperm.getMapping(j+1);
                                            int newj = Math.abs(map)-1;
                                            int newh = (map > 0 ? h : 1-h);
                                            MHNThrow sst2 = th[newj][newh][k+sss.getDelay()][slot];
                                            if (sst2 == null)
                                                throw new JuggleExceptionUser(errorstrings.getString("Error_badpattern_paths"));
                                            if (sst.pathnum == 0 || sst2.pathnum == 0)
                                                throw new JuggleExceptionUser(errorstrings.getString("Error_badpattern_paths"));
                                            if (pathmap[sst.pathnum] == 0)
                                                pathmap[sst.pathnum] = sst2.pathnum;
                                            else if (pathmap[sst.pathnum] != sst2.pathnum)
                                                throw new JuggleExceptionUser(errorstrings.getString("Error_badpattern_switchdelay"));
                                        }
                                    }
                                }
                            }
                        }
                    }
                        break;
                default:
                    throw new JuggleExceptionUser(errorstrings.getString("Error_unknown_symmetry"));
            }
            // convert path mapping to a string
            String pathmapstring = "";
            for (int j = 1; j < balls; j++)
                pathmapstring += pathmap[j] + ",";
            pathmapstring += pathmap[balls];

            JMLSymmetry sym = new JMLSymmetry(symtype, sss.getNumberOfJugglers(),
                                                sss.getJugglerPerm().toString(),
                                                getNumberOfPaths(), pathmapstring,
                                                (double)sss.getDelay() / bps);

            result.addSymmetry(sym);
        }

        /*
        Permutation delayperm = null;
        for (int z = 0; z < result.getNumberOfSymmetries(); z++) {      // store delay permutation for later
            JMLSymmetry tempsym = result.getSymmetry(z);
            if (tempsym.getType() == JMLSymmetry.TYPE_DELAY)
                delayperm = tempsym.getPathPerm();
        }*/

        // Step 3 -- Add the primary events to the pattern:

        // We'll need to keep track of which hands/paths don't get any events,
        // so we can add positioning events in later steps
        // boolean[] pathcaught = new boolean[p.getNumberOfPaths()];
        boolean[][] handtouched = new boolean[getNumberOfJugglers()][2];
        for (int j = 0; j < getNumberOfJugglers(); j++)
            for (int h = 0; h < 2; h++)
                handtouched[j][h] = false;
        boolean[] pathtouched = new boolean[getNumberOfPaths()];
        for (int j = 0; j < getNumberOfPaths(); j++)
            pathtouched[j] = false;

        MHNThrow[][][][] th = getThrows();

        for (int k = 0; k < getPeriod(); k++) {
            for (int j = 0; j < getNumberOfJugglers(); j++) {
                for (int h = 0; h < 2; h++) {

                    MHNThrow sst = th[j][h][k][0];
                    if (sst == null || sst.master != sst)
                        continue;

                    // Step 3a -- Add transitions to the on-beat event (throw or holding transitions):

                    JMLEvent ev = new JMLEvent();
                    double throwxsum = 0.0;
                    int num_throws = 0;
                    boolean onethrown = false;

                    for (int slot = 0; slot < getMaxOccupancy(); slot++) {
                        MHNThrow sst2 = th[j][h][k][slot];
                        if (sst2 == null)
                            continue;

                        String type = null, mod = null;
                        switch (sst2.mod.charAt(0)) {
                            case 'B':
                                type = "bounce";
                                mod = null;
                                if (sst2.mod.indexOf("F") != -1) {
                                    if (mod == null)
                                        mod = "forced=true";
                                    else
                                        mod = mod + ";forced=true";
                                }
                                if (sst2.mod.indexOf("H") != -1) {
                                    if (mod == null)
                                        mod = "hyper=true";
                                    else
                                        mod = mod + ";hyper=true";
                                }
                                int bounces = 1;
                                for (int i = 1; i < sst2.mod.length(); i++) {
                                    if (sst2.mod.charAt(i) == 'B')
                                        bounces++;
                                }
                                if (bounces > 1) {
                                    if (mod == null)
                                        mod = "bounces=" + bounces;
                                    else
                                        mod = mod + ";bounces=" + bounces;
                                }
                                if (bouncefrac != MHNPattern.bouncefrac_default)
                                    if (mod == null)
                                        mod = "bouncefrac=" + bouncefrac;
                                    else
                                        mod = mod + ";bouncefrac=" + bouncefrac;
                                if (gravity != MHNPattern.gravity_default) {
                                    if (mod == null)
                                        mod = "g=" + gravity;
                                    else
                                        mod = mod + ";g=" + gravity;
                                }
                                break;
                            case 'F':
                                type = "bounce";
                                mod = "forced=true";
                                if (bouncefrac != MHNPattern.bouncefrac_default)
                                    mod += ";bouncefrac="+bouncefrac;
                                if (gravity != MHNPattern.gravity_default)
                                    mod = mod + ";g=" + gravity;
                                break;
                            case 'H':
                                type = "hold";
                                mod = null;
                                break;
                            case 'T':
                            default:
                                type = "toss";
                                mod = null;
                                if (gravity != MHNPattern.gravity_default)
                                    mod = "g=" + gravity;
                                break;
                        }

                        if (sst2.mod.charAt(0) != 'H') {
                            ev.addTransition(new JMLTransition(JMLTransition.TRANS_THROW, sst2.pathnum, type, mod));

                            int throwval = sst2.targetindex - k;
                            // System.out.println("slot = "+k+", hand = "+h+", throwval = "+throwval);

                            if (throwval == 1)
                                onethrown = true;
                            if (sst2.targethand == h) {
                                throwxsum += (throwval > 8 ? samethrowx[8] : samethrowx[throwval]);
                            } else {
                                throwxsum += (throwval > 8 ? crossingthrowx[8] : crossingthrowx[throwval]);
                            }
                            num_throws++;
                        } else if (hands != null) {
                            if (sst2.pathnum != -1) {   // -1 signals a 0 throw
                                // add holding transition if there's a ball in hand and "hands" is specified
                                ev.addTransition(new JMLTransition(JMLTransition.TRANS_HOLDING, sst2.pathnum, type, mod));
                                pathtouched[sst2.pathnum-1] = true;
                            }
                        }
                    }

                    // Step 3b -- Finish off the on-beat event based on the transitions we've added:

                    // set the event position
                    if (hands == null) {
                        if (num_throws > 0) {
                            double throwxav = throwxsum / (double)num_throws;
                            if (h == MHNPattern.LEFT_HAND)
                                throwxav = -throwxav;
                            ev.setLocalCoordinate(new Coordinate(throwxav,0.0,0.0));
                            ev.calcpos = false;
                        } else {
                            // mark event to calculate coordinate later
                            ev.calcpos = true;
                        }
                    } else {
                        Coordinate c = hands.getCoordinate(sst.juggler, sst.handsindex, 0);
                        if (h == MHNPattern.LEFT_HAND)
                            c.x = -c.x;
                        ev.setLocalCoordinate(c);
                        ev.calcpos = false;
                    }

                    // set the event time
                    double throwtime;
                    if (onethrown)
                        throwtime = ((double)k - 0.25*dwell) / bps;
                    else
                        throwtime = (double)k / bps;
                    ev.setT(throwtime);

                    // set the event juggler and hand
                    ev.setHand(j+1, (h==MHNPattern.RIGHT_HAND ? HandLink.RIGHT_HAND : HandLink.LEFT_HAND));

                    // add it to the pattern
                    result.addEvent(ev);

                    // record which hands are touched by this event, for later reference
                    for (int i2 = 0; i2 < indexes; i2++) {
                        for (int j2 = 0; j2 < numjugglers; j2++) {
                            for (int h2 = 0; h2 < 2; h2++) {
                                for (int slot2 = 0; slot2 < max_occupancy; slot2++) {
                                    MHNThrow sst2 = th[j2][h2][i2][slot2];
                                    if (sst2 != null && sst2.master == sst)
                                        handtouched[j2][h2] = true;
                                }
                            }
                        }
                    }


                    // Step 3c -- Add any catching (or holding) events immediately prior to the on-beat event
                    // added in Step 3b above:

                    // calculate pathcaught[], catchxsum, num_catches, and onecaught
                    /* for (int z = 0; z < p.getNumberOfPaths(); z++)
                        pathcaught[z] = false; */
                    double catchxsum = 0.0;
                    int num_catches = 0;
                    boolean onecaught = false;

                    for (int slot = 0; slot < getMaxOccupancy(); slot++) {
                        MHNThrow sst2 = th[j][h][k][slot];
                        if (sst2 == null)
                            break;
                        if (!sst2.catching)
                            continue;

                        int catchpath = sst2.pathnum;
                        int catchval = k - sst2.source.index;
                        /* if (pathcaught[catchpath-1] == true)
                            throw new JuggleExceptionInternal("Caught path "+catchpath+" twice");*/
                        // pathcaught[catchpath-1] = true;
                        pathtouched[catchpath-1] = true;
                        // System.out.println("catching value "+catchval);
                        catchxsum += (catchval > 8 ? catchx[8] : catchx[catchval]);
                        num_catches++;
                        if (catchval == 1)
                            onecaught = true;
                    }

                    // Now add the event(s).  There are two cases to consider:  (1) all catches happen at
                    // the same event, or (2) multiple catch events are made in succession.
                    double lastcatchtime = 0.0;

                    if (splitcatchfactor == 0.0 || num_catches < 2) {
                        // Case 1: everything happens at a single event
                        ev = new JMLEvent();

                        // first set the event position
                        if (hands == null) {
                            if (num_catches > 0) {
                                double cx = catchxsum / (double)num_catches;
                                // System.out.println("average catch pos. = "+cx);
                                ev.setLocalCoordinate(new Coordinate((h==MHNPattern.RIGHT_HAND?cx:-cx),0.0,0.0));
                                ev.calcpos = false;
                            } else {
                                // mark event to calculate coordinate later
                                ev.calcpos = true;
                            }
                        } else {
                            int pos = sst.handsindex - 2;
                            while (pos < 0)
                                pos += hands.getPeriod(sst.juggler);
                            int index = hands.getCatchIndex(sst.juggler, pos);
                            Coordinate c = hands.getCoordinate(sst.juggler, pos, index);
                            if (h == MHNPattern.LEFT_HAND)
                                c.x = -c.x;
                            ev.setLocalCoordinate(c);
                            ev.calcpos = false;
                        }

                        // set the event time
                        if (onecaught)
                            lastcatchtime = ((double)k - 0.5*dwell) / bps;
                        else
                            lastcatchtime = ((double)k - dwell) / bps;

//Mahit begin
                        if (hss != null && dwellmax) {
                        	lastcatchtime = ((double)k - dwellarray[k]) / bps;
                        }
//Mahit end                        
                        ev.setT(lastcatchtime);

                        // set the event juggler and hand
                        ev.setHand(j+1, (h==MHNPattern.RIGHT_HAND ? HandLink.RIGHT_HAND : HandLink.LEFT_HAND));

                        // add all the transitions
                        for (int slot = 0; slot < getMaxOccupancy(); slot++) {
                            MHNThrow sst2 = th[j][h][k][slot];
                            if (sst2 == null)
                                break;

                            if (sst2.catching) {
                                ev.addTransition(new JMLTransition(JMLTransition.TRANS_CATCH, sst2.pathnum, null, null));
                            } else if (hands != null) {
                                if (sst2.pathnum != -1) {   // -1 signals a 0 throw
                                    // add holding transition if there's a ball in hand and "hands" is specified
                                    ev.addTransition(new JMLTransition(JMLTransition.TRANS_HOLDING, sst2.pathnum, null, null));
                                    pathtouched[sst2.pathnum-1] = true;
                                }
                            }
                        }

                        // add event to the pattern
                        result.addEvent(ev);
                    } else {
                        // Case 2: separate event for each catch; we know that numcatches>1 here
                        for (int slot = 0; slot < getMaxOccupancy(); slot++) {
                            MHNThrow sst2 = th[j][h][k][slot];
                            if (sst2 == null)
                                continue;

                            // we only need to add catch transitions here; holding transitions will be
                            // added in Step 8
                            if (sst2.catching) {
                                ev = new JMLEvent();

                                // first set the event position
                                if (hands == null) {
                                    double cx = catchxsum / (double)num_catches;
                                    // System.out.println("average catch pos. = "+cx);
                                    ev.setLocalCoordinate(new Coordinate((h==MHNPattern.RIGHT_HAND?cx:-cx),0.0,0.0));
                                } else {
                                    int pos = sst.handsindex - 2;
                                    while (pos < 0)
                                        pos += hands.getPeriod(sst.juggler);
                                    int index = hands.getCatchIndex(sst.juggler, pos);
                                    Coordinate c = hands.getCoordinate(sst.juggler, pos, index);
                                    if (h == MHNPattern.LEFT_HAND)
                                        c.x = -c.x;
                                    ev.setLocalCoordinate(c);
                                }
                                ev.calcpos = false;

                                // set the event time
                                double catchtime;
                                if (onecaught)
                                    catchtime = ((double)k - 0.5*dwell + ((double)sst2.catchnum/(double)(num_catches-1) - 0.5) *
                                                 splitcatchfactor) / bps;
                                else
                                    catchtime = ((double)k - dwell + ((double)sst2.catchnum/(double)(num_catches-1) - 0.5) *
                                                 splitcatchfactor) / bps;
                                ev.setT(catchtime);

                                if (sst.catchnum == (num_catches - 1))
                                    lastcatchtime = catchtime;

                                // set the event juggler and hand
                                ev.setHand(j+1, (h==MHNPattern.RIGHT_HAND ? HandLink.RIGHT_HAND : HandLink.LEFT_HAND));

                                // add the transition
                                ev.addTransition(new JMLTransition(JMLTransition.TRANS_CATCH, sst2.pathnum, null, null));

                                // add catch event to the pattern
                                result.addEvent(ev);
                            }
                        }
                    }


                    // Step 3d -- Add other hand positioning events between the catch/throw/catch events,
                    // if they are specified:

                    if (hands == null)
                        continue;

                    // add other events between the previous catch and the current throw
                    int pos = sst.handsindex - 2;
                    while (pos < 0)
                        pos += hands.getPeriod(sst.juggler);
                    int catchindex = hands.getCatchIndex(sst.juggler, pos);
                    int numcoords = hands.getNumberOfCoordinates(sst.juggler, pos) - catchindex;

                    for (int di = 1; di < numcoords; di++) {
                        Coordinate c = hands.getCoordinate(sst.juggler, pos, catchindex+di);
                        if (c == null)
                            continue;
                        ev = new JMLEvent();
                        if (h == MHNPattern.LEFT_HAND)
                            c.x = -c.x;
                        ev.setLocalCoordinate(c);
                        ev.calcpos = false;
                        ev.setT(lastcatchtime + (double)di*(throwtime-lastcatchtime)/numcoords);
                        ev.setHand(sst.juggler, (h==MHNPattern.RIGHT_HAND?HandLink.RIGHT_HAND:HandLink.LEFT_HAND));
                        result.addEvent(ev);
                    }

                    // figure out when the next catch or hold is:
                    double nextcatchtime = lastcatchtime;

                    for (int tempk = (k+1); tempk < getIndexes(); tempk++) {
                        int next_num_catches = 0;
                        boolean next_gotevent = false;
                        boolean next_onecaught = false;

                        for (int tempslot = 0; tempslot < getMaxOccupancy(); tempslot++) {
                            MHNThrow tempsst = th[j][h][tempk][tempslot];
                            if (tempsst == null)
                                break;

                            next_gotevent = true;

                            if (tempsst.catching) {
                                next_num_catches++;
                                if ((tempk - tempsst.source.index) == 1)
                                    next_onecaught = true;
                            }
                        }

                        if (next_gotevent) {
                            if (next_num_catches < 2) {
                                if (next_onecaught)
                                    nextcatchtime = ((double)tempk - 0.5*dwell) / bps;
                                else
                                    nextcatchtime = ((double)tempk - dwell) / bps;
                            } else {
                                if (next_onecaught)
                                    nextcatchtime = ((double)tempk - 0.5*dwell - 0.5*splitcatchfactor) / bps;
                                else
                                    nextcatchtime = ((double)tempk - dwell - 0.5*splitcatchfactor) / bps;
                            }
                            break;
                        }
                    }
                    if (nextcatchtime == lastcatchtime)
                        throw new JuggleExceptionInternal("Couldn't find next catch/hold past t="+lastcatchtime);

                    // add other events between the current throw and the next catch
                    pos = sst.handsindex;
                    numcoords = hands.getCatchIndex(sst.juggler, pos);

                    for (int di = 1; di < numcoords; di++) {
                        Coordinate c = hands.getCoordinate(sst.juggler, pos, di);
                        if (c == null)
                            continue;
                        ev = new JMLEvent();
                        if (h == MHNPattern.LEFT_HAND)
                            c.x = -c.x;
                        ev.setLocalCoordinate(c);
                        ev.calcpos = false;
                        ev.setT(throwtime + (double)di*(nextcatchtime-throwtime)/numcoords);
                        ev.setHand(sst.juggler, (h==MHNPattern.RIGHT_HAND?HandLink.RIGHT_HAND:HandLink.LEFT_HAND));
                        result.addEvent(ev);
                    }
                }

                // Step 3e -- Define a body position for this juggler and beat, if one is specified:

                if (bodies != null) {
                    int index = k % bodies.getPeriod(j+1);
                    int coords = bodies.getNumberOfPositions(j+1, index);
                    for (int z = 0; z < coords; z++) {
                        JMLPosition jmlp = bodies.getPosition(j+1, index, z);
                        if (jmlp != null) {
                            jmlp.setT(((double)k + (double)z / (double)coords) / bps);
                            result.addPosition(jmlp);
                        }
                    }
                }
            }
        }

        // Step 4 -- Add simple positioning events for hands that got no events:

        for (int j = 0; j < getNumberOfJugglers(); j++) {
            for (int h = 0; h < 2; h++) {
                if (!handtouched[j][h]) {
                    JMLEvent ev = new JMLEvent();
                    ev.setLocalCoordinate(new Coordinate((h==MHNPattern.RIGHT_HAND?restingx:-restingx),0.0,0.0));
                    ev.setT(-1.0);
                    ev.setHand(j+1, (h==0?HandLink.RIGHT_HAND:HandLink.LEFT_HAND));
                    ev.calcpos = false;
                    result.addEvent(ev);
                }
            }
        }

        // Step 5 -- Add <holding> transitions for paths that got no events:

        //     first, apply all pattern symmetries to figure out which paths don't get touched
        for (int j = 0; j < result.getNumberOfSymmetries(); j++) {
            Permutation perm = result.getSymmetry(j).getPathPerm();
            for (int k = 0; k < getNumberOfPaths(); k++) {
                if (pathtouched[k]) {
                    for (int l = 1; l < perm.getOrder(k+1); l++) {
                        pathtouched[perm.getMapping(k+1,l) - 1] = true;
                    }
                }
            }
        }
        //     next, add <holding> transitions for the untouched paths
        for (int k = 0; k < getNumberOfPaths(); k++) {
            if (pathtouched[k])
                continue;

            // figure out which hand it should belong in
            int hand = HandLink.LEFT_HAND;
            int juggler = 0;

top:
            for (int tempk = 0; tempk < getIndexes(); tempk++) {
                for (int tempj = 0; tempj < getNumberOfJugglers(); tempj++) {
                    for (int temph = 0; temph < 2; temph++) {
                        for (int slot = 0; slot < getMaxOccupancy(); slot++) {
                            MHNThrow sst = th[tempj][temph][tempk][slot];
                            if (sst != null && sst.pathnum == (k+1)) {
                                hand = (temph==MHNPattern.RIGHT_HAND?HandLink.RIGHT_HAND:HandLink.LEFT_HAND);
                                juggler = tempj;
                                break top;
                            }
                        }
                    }
                }
            }

            // add <holding> transitions to each of that hand's events
            JMLEvent ev = result.getEventList();
            while (ev != null) {
                if (ev.getHand() == hand && ev.getJuggler() == (juggler+1)) {
                    ev.addTransition(new JMLTransition(JMLTransition.TRANS_HOLDING, (k+1), null, null));

                    // mark related paths as touched
                    pathtouched[k] = true;
                    for (int j = 0; j < result.getNumberOfSymmetries(); j++) {
                        Permutation perm = result.getSymmetry(j).getPathPerm();
                        for (int l = 1; l < perm.getOrder(k+1); l++) {
                            pathtouched[perm.getMapping(k+1,l) - 1] = true;
                        }
                    }
                }

                ev = ev.getNext();
            }
            //  if (ev == null)
            //      throw new JuggleExceptionUser("Could not find event for hand");
        }

        // Step 6 -- Do a build of the full event list so we can scan through it chronologically in Steps 7 and 8:

        result.buildEventList();

        // Step 7 -- Specify positions for events that don't have them defined yet:

        for (int j = 1; j <= getNumberOfJugglers(); j++) {
            for (int h = 0; h < 2; h++) {
                int hand = (h==MHNPattern.RIGHT_HAND?HandLink.RIGHT_HAND:HandLink.LEFT_HAND);

                JMLEvent ev = result.getEventList();
                JMLEvent start = null;
                int scanstate = 1;  // 1 = starting, 2 = on defined event, 3 = on undefined event
                while (ev != null) {
                    if (ev.getJuggler() == j && ev.getHand() == hand) {
                        // System.out.println("j = "+j+", h = "+h+", t = "+ev.getT()+", calcpos = "+ev.calcpos);
                        // System.out.println("   initial state = "+scanstate);

                        if (ev.calcpos) {
                            scanstate = 3;
                        } else {
                            switch (scanstate) {
                                case 1:
                                    scanstate = 2;
                                    break;
                                case 2:
                                    break;
                                case 3:
                                    if (start != null) {
                                        JMLEvent end = ev;
                                        double t_end = end.getT();
                                        Coordinate pos_end = end.getLocalCoordinate();
                                        double t_start = start.getT();
                                        Coordinate pos_start = start.getLocalCoordinate();

                                        ev = start.getNext();
                                        while (ev != end) {
                                            if (ev.getJuggler() == j && ev.getHand() == hand) {
                                                double t = ev.getT();
                                                double x = pos_start.x + (t - t_start)*
                                                    (pos_end.x - pos_start.x) / (t_end - t_start);
                                                double y = pos_start.y + (t - t_start)*
                                                    (pos_end.y - pos_start.y) / (t_end - t_start);
                                                double z = pos_start.z + (t - t_start)*
                                                    (pos_end.z - pos_start.z) / (t_end - t_start);
                                                ev.setLocalCoordinate(new Coordinate(x, y, z));
                                                ev.calcpos = false;
                                            }
                                            ev = ev.getNext();
                                        }
                                    }
                                    scanstate = 2;
                                    break;
                            }
                            start = ev;
                        }
                        // System.out.println("   final state = "+scanstate);
                    }
                    ev = ev.getNext();
                }

                // do a last scan through to define any remaining undefined positions
                ev = result.getEventList();
                while (ev != null) {
                    if (ev.getJuggler() == j && ev.getHand() == hand && ev.calcpos) {
                        ev.setLocalCoordinate(new Coordinate((h==MHNPattern.RIGHT_HAND?restingx:-restingx),0.0,0.0));
                        ev.calcpos = false;
                    }
                    ev = ev.getNext();
                }
            }
        }

        // Step 8 -- Scan through the list of events, and look for cases where we need to add additional
        // <holding> transitions.  These are marked by cases where the catch and throw transitions for
        // a given path have intervening events in that hand; we want to add <holding> transitions to
        // these intervening events:

        for (int k = 0; k < getNumberOfPaths(); k++) {
            boolean add_mode = false;
            boolean found_event = false;
            int add_juggler = 0, add_hand = 0;

            JMLEvent ev = result.getEventList();
            while (ev != null) {
                JMLTransition tr = ev.getPathTransition((k+1), JMLTransition.TRANS_ANY);
                if (tr != null) {
                    switch (tr.getType()) {
                        case JMLTransition.TRANS_THROW:
                            if (!found_event && !add_mode) {
                                // first event mentioning path is a throw
                                // rewind to beginning of list and add holds
                                add_mode = true;
                                add_juggler = ev.getJuggler();
                                add_hand = ev.getHand();
                                ev = result.getEventList();
                                continue;
                            }
                            add_mode = false;
                            break;
                        case JMLTransition.TRANS_CATCH:
                        case JMLTransition.TRANS_SOFTCATCH:
                            add_mode = true;
                            add_juggler = ev.getJuggler();
                            add_hand = ev.getHand();
                            break;
                        case JMLTransition.TRANS_HOLDING:
                            if (!found_event && !add_mode) {
                                // first event mentioning path is a hold
                                // rewind to beginning of list and add holds
                                add_mode = true;
                                add_juggler = ev.getJuggler();
                                add_hand = ev.getHand();
                                ev = result.getEventList();
                                continue;
                            }
                            add_mode = true;
                            add_juggler = ev.getJuggler();
                            add_hand = ev.getHand();
                            break;
                    }
                    found_event = true;
                } else if (add_mode) {
                    if (ev.getJuggler()==add_juggler && ev.getHand()==add_hand && ev.isMaster())
                        ev.addTransition(new JMLTransition(JMLTransition.TRANS_HOLDING, (k+1), null, null));
                }
                ev = ev.getNext();
            }
        }

        if (title != null)
            result.setTitle(title);

        if (Constants.DEBUG_LAYOUT) {
            System.out.println("Pattern in JML format:\n");
            System.out.println(result);
        }

        return result;
    }

    protected static final double[] throwspersec =
        { 2.00, 2.00, 2.00, 2.90, 3.40, 4.10, 4.25, 5.00, 5.00, 5.50 };

    protected double calcBps() {
        // Calculate a default beats per second (bps) for the pattern
        double result = 0.0;

        MHNThrow[][][][] th = getThrows();
        int numberaveraged = 0;

        for (int k = 0; k < getPeriod(); k++) {
            for (int j = 0; j < getNumberOfJugglers(); j++) {
                for (int h = 0; h < 2; h++) {
                    for (int slot = 0; slot < getMaxOccupancy(); slot++) {
                        MHNThrow sst = th[j][h][k][slot];
                        if (sst != null) {
                            int throwval = sst.targetindex - k;
                            if (throwval > 2) {
                                result += throwspersec[throwval > 9 ? 9 : throwval];
                                numberaveraged++;
                            }
                        }
                    }
                }
            }
        }
        if (numberaveraged > 0)
            result /= (double)numberaveraged;
        else
            result = 2.0;

        return result;
    }

}
