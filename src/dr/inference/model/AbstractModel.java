/*
 * AbstractModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.model;

import dr.inference.parallel.MPISerializable;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.ArrayList;
import java.util.Map;

/**
 * A model that brings together a number of model components
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: AbstractModel.java,v 1.13 2006/08/17 15:30:08 rambaut Exp $
 */
public abstract class AbstractModel implements Model, ModelListener, VariableListener, StatisticList, MPISerializable {

    /**
     * @param name Model Name
     */
    public AbstractModel(String name) {
        this.name = name;
    }

    /**
     * Adds a sub-model to this model. If the model is already in the
     * list then it does nothing.
     */
    public void addModel(Model model) {
        Model.CONNECTED_MODEL_SET.add(model);

        if (!models.contains(model)) {
            models.add(model);
            model.addModelListener(this);
        }
    }

    public void removeModel(Model model) {
        models.remove(model);
        model.removeModelListener(this);
    }

    public int getModelCount() {
        return models.size();
    }

    public final Model getModel(int i) {
        return models.get(i);
    }

    public final void addParameter(Parameter parameter) {
        if (parameter instanceof Parameter) {
            Parameter.CONNECTED_PARAMETER_SET.add((Parameter) parameter);
        }

        if (!parameters.contains(parameter)) {
            parameters.add(parameter);
            parameter.addVariableListener(this);
        }

        // parameters are also statistics
        if (parameter instanceof Statistic) addStatistic((Statistic) parameter);
    }

    public final void removeParameter(Parameter parameter) {
        parameters.remove(parameter);
        parameter.removeVariableListener(this);

        // parameters are also statistics
        if (parameter instanceof Statistic) removeStatistic((Statistic) parameter);
    }

    /**
     * @param parameter
     * @return true of the given parameter is contained in this model
     */
    public final boolean hasParameter(Parameter parameter) {
        return parameters.contains(parameter);
    }

    /**
     * Adds a model listener.
     */
    public void addModelListener(ModelListener listener) {
        listenerHelper.addModelListener(listener);
    }

    /**
     * remove a model listener.
     */
    public void removeModelListener(ModelListener listener) {
        listenerHelper.removeModelListener(listener);
    }

    public void addModelRestoreListener(ModelListener listener) {
        listenerHelper.addModelRestoreListener(listener);
    }

    public boolean isUsed() {
        return listenerHelper.getListenerCount() > 0;
    }


    /**
     * Fires a model changed event.
     */
    public void fireModelChanged() {
        listenerHelper.fireModelChanged(this, this, -1);
    }

    public void fireModelChanged(Object object) {
        listenerHelper.fireModelChanged(this, object, -1);
    }

    public void fireModelChanged(Object object, int index) {
        listenerHelper.fireModelChanged(this, object, index);
    }

    public final int getParameterCount() {
        return parameters.size();
    }

    public final Variable getVariable(int i) {
        return parameters.get(i);
    }

    // **************************************************************
    // MPI IMPLEMENTATION
    // **************************************************************


    public void sendState(int toRank) {

        // Iterate through child models
        for (Model model : models) {
            ((AbstractModel) model).sendState(toRank);
        }
        // Send current model parameters
        for (Parameter parameter : parameters) {
            if (parameter instanceof Parameter.Abstract) ((Parameter.Abstract) parameter).sendState(toRank);
        }
    }


    public void sendStateNoParameters(int toRank) {

        // Iterate through child models
        for (Model model : models) {
            ((AbstractModel) model).sendState(toRank);
        }
    }

    public void receiveStateNoParameters(int fromRank) {
        for (Model model : models) {
            ((AbstractModel) model).receiveState(fromRank);
        }
    }

    public void receiveState(int fromRank) {
        for (Model model : models) {
            ((AbstractModel) model).receiveState(fromRank);
        }
        // Send current model parameters
        for (Parameter parameter : parameters) {
            if (parameter instanceof Parameter.Abstract)
                ((Parameter.Abstract) parameter).receiveState(fromRank);
        }


    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    public final void modelChangedEvent(Model model, Object object, int index) {

//		String message = "  model: " + getModelName() + "/" + getId() + "  component: " + model.getModelName();
//		if (object != null) {
//			message += " object: " + object;
//		}
//		if (index != -1) {
//			message += " index: " + index;
//		}
//		System.out.println(message);

        handleModelChangedEvent(model, object, index);
    }

    // do nothing by default
    public void modelRestored(Model model) {
    }

    abstract protected void handleModelChangedEvent(Model model, Object object, int index);

    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************

    public final void variableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        handleVariableChangedEvent(variable, index, type);

        // todo AR - I am not sure this is required and may be overruling modelChange events on parts of the
        // model. If a parameter changes it should be handleVariableChangedEvent() job to fireModelChanged
        // events
        listenerHelper.fireModelChanged(this, variable, index);
    }

    /**
     * This method is called whenever a parameter is changed.
     * <p/>
     * It is strongly recommended that the model component sets a "dirty" flag and does no
     * further calculations. Recalculation is typically done when the model component is asked for
     * some information that requires them. This mechanism is 'lazy' so that this method
     * can be safely called multiple times with minimal computational cost.
     */
    protected abstract void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type);

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    public final void storeModelState() {
        if (isValidState) {
            //System.out.println("STORE MODEL: " + getModelName() + "/" + getId());

            for (Model m : models) {
                m.storeModelState();
            }

            for (Parameter parameter : parameters) {
                parameter.storeModelState();
            }

            storeState();
            isValidState = false;
        }
    }

    public final void restoreModelState() {
        if (!isValidState) {
            //System.out.println("RESTORE MODEL: " + getModelName() + "/" + getId());

            for (Parameter parameter : parameters) {
                parameter.restoreModelState();
            }
            for (Model m : models) {
                m.restoreModelState();
            }

            restoreState();
            isValidState = true;

            listenerHelper.fireModelRestored(this);
        }
    }

    public final void acceptModelState() {
        if (!isValidState) {
            //System.out.println("ACCEPT MODEL: " + getModelName() + "/" + getId());

            for (Model m : models) {
                m.acceptModelState();
            }

            acceptState();
            isValidState = true;
        }
    }

    @Override
    public void saveModelState(Map<String, Object> stateMap) {
        if (!isValidState) {
            //System.out.println("ACCEPT MODEL: " + getModelName() + "/" + getId());

            for (Parameter parameter : parameters) {
                parameter.saveModelState(stateMap);
            }

            for (Model m : models) {
                m.saveModelState(stateMap);
            }

            saveState(stateMap);
        }
    }

    @Override
    public void loadModelState(Map<String, Object> stateMap) {
        if (!isValidState) {
            //System.out.println("ACCEPT MODEL: " + getModelName() + "/" + getId());

            for (Parameter parameter : parameters) {
                parameter.saveModelState(stateMap);
            }

            for (Model m : models) {
                m.saveModelState(stateMap);
            }

            loadState(stateMap);

            isValidState = true;
        }
    }

    public boolean isValidState() {
        return isValidState;
    }

    public final String getModelName() {
        return name;
    }


    /**
     * Additional state information, outside of the sub-model is stored by this call.
     */
    protected abstract void storeState();

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    protected abstract void restoreState();

    /**
     * After this call the model is guaranteed to have returned its extra state information to
     * the values coinciding with the last storeState call.
     * Sub-models are handled automatically and do not need to be considered in this method.
     */
    protected void acceptState() {
        // do nothing by default - override
    }

    protected void saveState(Map<String, Object> stateMap) {
        throw new UnsupportedOperationException("must be overridden");
    }

    protected void loadState(Map<String, Object> loadMap) {
        throw new UnsupportedOperationException("must be overridden");
    }

    // **************************************************************
    // StatisticList IMPLEMENTATION
    // **************************************************************

    public final void addStatistic(Statistic statistic) {
        if (!statistics.contains(statistic)) {
            statistics.add(statistic);
        }
    }

    public final void removeStatistic(Statistic statistic) {

        statistics.remove(statistic);
    }

    /**
     * @return the number of statistics of this component.
     */
    public int getStatisticCount() {

        return statistics.size();
    }

    /**
     * @return the ith statistic of the component
     */
    public Statistic getStatistic(int i) {

        return statistics.get(i);
    }

    public final Statistic getStatistic(String name) {

        for (int i = 0; i < getStatisticCount(); i++) {
            Statistic statistic = getStatistic(i);
            if (name.equals(statistic.getStatisticName())) {
                return statistic;
            }
        }

        return null;
    }

    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

    private String id = null;

    public void setId(String id) {
        this.id = id;
    }

    public String getId() {
        return id;
    }


    public String toString() {
        if (id != null) {
            return id;
        } else if (name != null) {
            return name;
        }
        return super.toString();
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document d) {
        throw new RuntimeException("Not implemented!");
    }

    boolean isValidState = true;

    protected Model.ListenerHelper listenerHelper = new Model.ListenerHelper();

    private final ArrayList<Model> models = new ArrayList<Model>();
    private final ArrayList<Parameter> parameters = new ArrayList<Parameter>();
    private final ArrayList<Variable> variables = new ArrayList<Variable>();
    private final ArrayList<Statistic> statistics = new ArrayList<Statistic>();

    private final String name;
}
