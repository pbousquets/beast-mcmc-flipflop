package dr.evomodelxml.substmodel;

import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.AFsequence;
import dr.evolution.datatype.DataType;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evoxml.util.DataTypeUtils;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.text.NumberFormat;
import java.util.logging.Logger;

public class FlipFlopCenancestorFrequencyParser extends AbstractXMLObjectParser {
    public static final String FREQUENCIES = "frequencies";
    public static final String FREQUENCY_MODEL = "flipflopCenancestorFrequency";
    public static final String NORMALIZE = "normalize";
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
        SubstitutionModel model = (SubstitutionModel) xo.getChild(SubstitutionModel.class);
        int nStates=model.getFrequencyModel().getFrequencyCount();

        Double methylatedProportion=null;
        if(xo.hasAttribute(METHYLATEDPROPORTION)){
            methylatedProportion=xo.getDoubleAttribute(METHYLATEDPROPORTION);
        }

        Parameter freqsParam = (Parameter) xo.getElementFirstChild(FREQUENCIES);

        StringBuilder sb = new StringBuilder("Creating cenancestor state frequencies model '" + freqsParam.getParameterName() + "': ");

        if(methylatedProportion == null & (freqsParam.getDimension() != nStates)){
            throw new XMLParseException("dimension of frequency parameter and number of sequence states don't match with "+METHYLATEDPROPORTION+" not specified");
        }

        // frequencies set programatically
        if(methylatedProportion != null) {
            sb.append("\n\tSetting frequencies programatically using the "+METHYLATEDPROPORTION+ " parameter: "+methylatedProportion.toString());
            sb.append("\n\tFrequencies ");

            if (methylatedProportion < 0 | methylatedProportion > 1) throw new XMLParseException(METHYLATEDPROPORTION + "must be between 0 and 1");
            freqsParam.setDimension(nStates);
            for (int i=0; i<freqsParam.getDimension(); i++){
                freqsParam.setParameterValueQuietly(i,0.0);
            }
            //DM TODO Are these positions always as expected? Should we implement methods to give us these positions?
            freqsParam.setParameterValueQuietly(nStates-1,methylatedProportion);
            freqsParam.setParameterValue(0,1-methylatedProportion);
        } else {sb.append("\n\tInitial frequencies ");}
        sb.append("= {");

        double sum = 0;
        for (int j = 0; j < freqsParam.getDimension(); j++) {
            sum += freqsParam.getParameterValue(j);
        }

        if (xo.getAttribute(NORMALIZE, false)) {
            for (int j = 0; j < freqsParam.getDimension(); j++) {
                if (sum != 0)
                    freqsParam.setParameterValue(j, freqsParam.getParameterValue(j) / sum);
                else
                    freqsParam.setParameterValue(j, 1.0 / freqsParam.getDimension());
            }
            sum = 1.0;
        }

        if (Math.abs(sum - 1.0) > 1e-8) {
            throw new XMLParseException("Frequencies do not sum to 1 (they sum to " + sum + ")");
        }


        NumberFormat format = NumberFormat.getNumberInstance();
        format.setMaximumFractionDigits(5);

        sb.append(format.format(freqsParam.getParameterValue(0)));
        for (int j = 1; j < freqsParam.getDimension(); j++) {
            sb.append(", ");
            sb.append(format.format(freqsParam.getParameterValue(j)));
        }
        sb.append("}");
        Logger.getLogger("dr.evomodel").info(sb.toString());

        return new FrequencyModel(dataType, freqsParam);
    }

    public String getParserDescription() {
        return "A model of cenancestor frequencies.";
    }

    public Class getReturnType() {
        return FrequencyModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PatternList.class, "Initial value", 0, 1), //to set the datatype
            new ElementRule(SubstitutionModel.class), //to get the number of states
            new ElementRule(FREQUENCIES, new XMLSyntaxRule[]{new ElementRule(Parameter.class)}), //frequency parameter
            new ElementRule(METHYLATEDPROPORTION,double.class,"Proportion of fully-methylated cells",true), //to re-initialize the frequencies programatically if desired
            AttributeRule.newBooleanRule(NORMALIZE, true), // to normalize specified frequencies
    };
}
