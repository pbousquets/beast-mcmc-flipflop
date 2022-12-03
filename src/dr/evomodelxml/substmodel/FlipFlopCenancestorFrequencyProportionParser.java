package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.FlipFlopCenancestorFrequencyProportionModel;
import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

public class FlipFlopCenancestorFrequencyProportionParser extends AbstractXMLObjectParser {
    public static final String FREQUENCIES = "frequencies";
    public static final String FREQUENCY_MODEL = "flipflopCenancestorFrequencyProportion";
    public static final String METHYLATEDPROPORTION = "methylatedProportion";

    public String[] getParserNames() {
        return new String[]{
                getParserName(), "beast_" + getParserName()
        };
    }

    public String getParserName() {
        return FREQUENCY_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        DataType dataType = ((PatternList) xo.getChild(PatternList.class)).getDataType();
        Parameter methylatedProportionParam= (Parameter) xo.getElementFirstChild(METHYLATEDPROPORTION);

        StringBuilder sb = new StringBuilder("Creating cenancestor state using a methylation-proportion model:");
        double methylatedProportion=methylatedProportionParam.getParameterValue(0);
        if (methylatedProportion <= 0 | methylatedProportion >= 1) throw new XMLParseException(METHYLATEDPROPORTION + "must be between 0 and 1");

        sb.append("\n\tMethylation-proportion parameter \""+methylatedProportionParam.getParameterName()+"\" = "+methylatedProportionParam.toString());

        Logger.getLogger("dr.evomodel").info(sb.toString());

        return new FlipFlopCenancestorFrequencyProportionModel(dataType, methylatedProportionParam);
    }

    public String getParserDescription() {
        return "A model of cenancestor frequencies determined by the proportion of fully-methylated or fully-demethylated cells";
    }

    public Class getReturnType() {
        return FlipFlopCenancestorFrequencyProportionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PatternList.class, "Initial value", 0, 1), //to set the datatype
            new ElementRule(METHYLATEDPROPORTION,new XMLSyntaxRule[]{new ElementRule(Parameter.class)},"Proportion of fully-methylated cells"), //Methylation proportion parameter
    };
}
