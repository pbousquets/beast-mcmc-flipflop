package dr.inference.model;

import java.util.Map;

/**
 * Created by max on 11/30/15.
 */
public class ElementWiseMatrixMultiplicationParameter extends MatrixParameter {
    private MatrixParameter[] paramList;

    public ElementWiseMatrixMultiplicationParameter(String name) {
        super(name);
    }

    public ElementWiseMatrixMultiplicationParameter(String name, MatrixParameter[] matList) {
        super(name);
        this.paramList =matList;
        for (MatrixParameter mat : matList) {
            mat.addVariableListener(this);
        }
    }

    @Override
    public double getParameterValue(int dim) {
        double prod=1;
        for (int i = 0; i < paramList.length ; i++) {
            prod=prod* paramList[i].getParameterValue(dim);
        }
        return prod;
    }

    public double getParameterValue(int row, int col){
        double prod=1;
        for (int i = 0; i < paramList.length ; i++) {
            prod=prod* paramList[i].getParameterValue(row,col);
        }
        return prod;
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, ChangeType type) {
        fireParameterChangedEvent(index, type);
    }

    @Override
    public int getDimension() {
        return paramList[0].getDimension();
    }

    @Override
    public int getColumnDimension() {
        return paramList[0].getColumnDimension();
    }

    public int getRowDimension(){
        return paramList[0].getRowDimension();
    }
}
