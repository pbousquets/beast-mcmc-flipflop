package dr.evomodel.epidemiology.casetocase;

import dr.evolution.tree.TreeAttributeProvider;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodelxml.tree.TreeLoggerParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

import java.util.ArrayList;
import java.util.List;


/**
 * @author Matthew Hall
 */
public class PartitionedTreeLoggerParser extends TreeLoggerParser {

    private XMLSyntaxRule[] rules;

    public PartitionedTreeLoggerParser(){
        rules = new XMLSyntaxRule[super.getSyntaxRules().length + 1];
        System.arraycopy(super.getSyntaxRules(), 0, rules, 0, rules.length-1);
        rules[rules.length-1] = new ElementRule(CaseToCaseTreeLikelihood.class);
    }
    public String getParserName() {
        return "logPartitionedTree";
    }

    public String getParserDescription() {
        return "Logs a partitioned tree (phylogenetic tree and transmission tree)";
    }

    public Class getReturnType() {
        return PartitionedTreeLogger.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return this.rules;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        List<TreeAttributeProvider> taps = new ArrayList<TreeAttributeProvider>();
        List<TreeTraitProvider> ttps = new ArrayList<TreeTraitProvider>();
                      /*
        parseXMLParameters(xo, taps, ttps);
        CaseToCaseTreeLikelihood c2cTL = (CaseToCaseTreeLikelihood)xo.getChild(CaseToCaseTreeLikelihood.class);

        TreeAttributeProvider[] treeAttributeProviders = new TreeAttributeProvider[taps.size()];
        taps.toArray(treeAttributeProviders);
        TreeTraitProvider[] treeTraitProviders = new TreeTraitProvider[ttps.size()];
        ttps.toArray(treeTraitProviders);

        PartitionedTreeLogger logger = new PartitionedTreeLogger(c2cTL, tree, branchRates,
                treeAttributeProviders, treeTraitProviders,
                formatter, logEvery, nexusFormat, sortTranslationTable, mapNames, format, condition);

        if (title != null) {
            logger.setTitle(title);
        }

        return logger;  */
        return null;
    }
}

