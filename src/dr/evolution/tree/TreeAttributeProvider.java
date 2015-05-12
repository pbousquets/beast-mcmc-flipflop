package dr.evolution.tree;

import java.io.Serializable;

/**
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public interface TreeAttributeProvider extends Serializable {

	String[] getTreeAttributeLabel();

	String[] getAttributeForTree(Tree tree);
}
