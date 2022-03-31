// DNoteNotationControl.java
//
// Copyright 2002-2022 Jack Boyce and the Juggling Lab contributors

package jugglinglab.notation;


public class DNoteNotationControl extends MHNNotationControl {

	private static final long serialVersionUID = -2342851836323967518L;

	@Override
    public Pattern newPattern() {
        return new DNotePattern();
    }
}
