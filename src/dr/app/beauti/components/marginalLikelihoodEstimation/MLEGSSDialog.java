/*
 * MLEGSSDialog.java
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

package dr.app.beauti.components.marginalLikelihoodEstimation;

import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.types.TreePriorType;
import dr.app.beauti.util.PanelUtils;
import dr.app.gui.components.WholeNumberField;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import java.awt.event.ActionEvent;
import java.awt.event.KeyEvent;


/**
 * @author Guy Baele
 */
public class MLEGSSDialog {

    private JFrame frame;

    private final OptionsPanel optionsPanel;

    private JLabel labelPathSteps, labelChainLength, labelLogEvery, labelLogFileName;
    private JLabel labelStepDistribution, labelTreeWorkingPrior, labelParameterWorkingPrior;

    private WholeNumberField pathStepsField = new WholeNumberField(1, Integer.MAX_VALUE);
    private WholeNumberField chainLengthField = new WholeNumberField(1, Integer.MAX_VALUE);
    private WholeNumberField logEveryField = new WholeNumberField(1, Integer.MAX_VALUE);

    private JTextArea logFileNameField = new JTextArea("MLE.log");

    JCheckBox operatorAnalysis = new JCheckBox("Print operator analysis");

    private JComboBox stepDistribution = new JComboBox();
    private JComboBox treeWorkingPrior = new JComboBox();
    //private JComboBox parameterWorkingPrior = new JComboBox();

    private MarginalLikelihoodEstimationOptions options;
    private BeautiOptions beautiOptions;

    private String description = "Settings for marginal likelihood estimation using GSS";

    public MLEGSSDialog(final JFrame frame, final MarginalLikelihoodEstimationOptions options, final BeautiOptions beautiOptions) {
        this.frame = frame;
        this.options = options;
        this.beautiOptions = beautiOptions;

        optionsPanel = new OptionsPanel(12, 12);

        optionsPanel.setOpaque(false);

        JTextArea mleInfo = new JTextArea("Set the options to perform marginal likelihood " +
                "estimation (MLE) using generalized stepping-stone sampling (GSS).");
        mleInfo.setColumns(56);
        PanelUtils.setupComponent(mleInfo);
        optionsPanel.addSpanningComponent(mleInfo);

        pathStepsField.setValue(100);
        pathStepsField.setColumns(16);
        pathStepsField.setMinimumSize(pathStepsField.getPreferredSize());
        labelPathSteps = optionsPanel.addComponentWithLabel("Number of stepping stones:", pathStepsField);
        pathStepsField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                options.pathSteps = pathStepsField.getValue();
            }
        });

        chainLengthField.setValue(1000000);
        chainLengthField.setColumns(16);
        chainLengthField.setMinimumSize(chainLengthField.getPreferredSize());
        labelChainLength = optionsPanel.addComponentWithLabel("Length of chains:", chainLengthField);
        chainLengthField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                options.mleChainLength = chainLengthField.getValue();
            }
        });

        optionsPanel.addSeparator();

        logEveryField.setValue(1000);
        logEveryField.setColumns(16);
        logEveryField.setMinimumSize(logEveryField.getPreferredSize());
        labelLogEvery = optionsPanel.addComponentWithLabel("Log likelihood every:", logEveryField);
        logEveryField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                options.mleLogEvery = logEveryField.getValue();
            }
        });

        optionsPanel.addSeparator();

        logFileNameField.setColumns(32);
        logFileNameField.setEditable(false);
        logFileNameField.setMinimumSize(logFileNameField.getPreferredSize());
        labelLogFileName = optionsPanel.addComponentWithLabel("Log file name:", logFileNameField);
        logFileNameField.addKeyListener(new java.awt.event.KeyListener() {
            public void keyTyped(KeyEvent e) {
            }

            public void keyPressed(KeyEvent e) {
            }

            public void keyReleased(KeyEvent e) {
                //options.mleFileName = logFileNameField.getText();
            }
        });

        optionsPanel.addSeparator();

        JTextArea betaInfo = new JTextArea("By default, the power posteriors are determined according to " +
                "evenly spaced quantiles of a Beta(0.3, 1.0) distribution, thereby estimating " +
                "more power posteriors close to the prior.");
        betaInfo.setColumns(56);
        PanelUtils.setupComponent(betaInfo);
        optionsPanel.addSpanningComponent(betaInfo);

        treeWorkingPrior.addItem("Product of exponential distributions");
        treeWorkingPrior.addItem("Matching coalescent model");
        treeWorkingPrior.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                String selection = (String)((JComboBox)e.getSource()).getSelectedItem();
                if (selection.equals("Matching coalescent model")) {
                    beautiOptions.logCoalescentEventsStatistic = false;
                } else {
                    beautiOptions.logCoalescentEventsStatistic = true;
                }
                TreePriorType treePrior = beautiOptions.getPartitionTreePriors().get(0).getNodeHeightPrior();
                boolean mcmAllowed = false;
                if (treePrior.equals(TreePriorType.CONSTANT) || treePrior.equals(TreePriorType.EXPONENTIAL)
                        || treePrior.equals(TreePriorType.EXPANSION) || treePrior.equals(TreePriorType.LOGISTIC)) {
                    mcmAllowed = true;
                }
                if (selection.equals("Matching coalescent model") && !mcmAllowed) {
                    JOptionPane.showMessageDialog(frame,
                            "The selected coalescent model can't be equipped with a \n" +
                                    "matching working prior. Reverting back to use a product \n" +
                                    "of exponential distributions.",
                            "Matching coalescent model warning",
                            JOptionPane.WARNING_MESSAGE);
                    treeWorkingPrior.setSelectedItem("Product of exponential distributions");
                    options.choiceTreeWorkingPrior = "Product of exponential distributions";
                } else {
                    options.choiceTreeWorkingPrior = selection;
                }
            }
        });
        labelTreeWorkingPrior = optionsPanel.addComponentWithLabel("Tree working prior", treeWorkingPrior);

        //parameterWorkingPrior.addItem("Normal KDE");
        //parameterWorkingPrior.addItem("Gamma KDE");
        //labelParameterWorkingPrior = optionsPanel.addComponentWithLabel("Parameter working prior", parameterWorkingPrior);

        stepDistribution.addItem("Beta");
        labelStepDistribution = optionsPanel.addComponentWithLabel("Stepping stone distribution:", stepDistribution);

        optionsPanel.addSeparator();

        optionsPanel.addComponent(operatorAnalysis);

        operatorAnalysis.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(ActionEvent e) {
                if (operatorAnalysis.isSelected()) {
                    options.printOperatorAnalysis = true;
                } else {
                    options.printOperatorAnalysis = false;
                }
            }
        });

        optionsPanel.addSeparator();

        JTextArea mleTutorial = new JTextArea("Additional information on marginal likelihood estimation in BEAST " +
                "can be found on http://beast.bio.ed.ac.uk/Model-selection");
        mleTutorial.setColumns(56);
        PanelUtils.setupComponent(mleTutorial);
        optionsPanel.addSpanningComponent(mleTutorial);

        JTextArea citationText = new JTextArea("Baele G, Lemey P, Suchard MA (2015) Genealogical working " +
                "distributions for Bayesian \nmodel testing with phylogenetic uncertainty [GSS Paper].");
        citationText.setColumns(45);
        optionsPanel.addComponentWithLabel("Citation:", citationText);

        optionsPanel.addSeparator();
        
        /*JPanel panel = new JPanel(new FlowLayout(FlowLayout.CENTER));
        panel.add(optionsPanel, BorderLayout.CENTER);
        panel.setOpaque(false);
        
        JScrollPane scrollPane = new JScrollPane(panel, JScrollPane.VERTICAL_SCROLLBAR_AS_NEEDED, JScrollPane.HORIZONTAL_SCROLLBAR_NEVER);
        scrollPane.setOpaque(false);
        scrollPane.setBorder(null);
        scrollPane.getViewport().setOpaque(false);

        add(scrollPane, BorderLayout.CENTER);*/

    }

    public int showDialog() {

        JOptionPane optionPane = new JOptionPane(optionsPanel,
                JOptionPane.QUESTION_MESSAGE,
                JOptionPane.OK_CANCEL_OPTION,
                null,
                null,
                null);
        optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        final JDialog dialog = optionPane.createDialog(frame, description);
        dialog.pack();

        dialog.setVisible(true);

        int result = JOptionPane.CANCEL_OPTION;
        Integer value = (Integer) optionPane.getValue();
        if (value != null && value != -1) {
            result = value;
        }

        return result;
    }

	/*public JComponent getExportableComponent() {
		return optionsPanel;
	}*/

    public void setFilenameStem(String fileNameStem, boolean addTxt) {
        logFileNameField.setText(fileNameStem + ".mle.log" + (addTxt ? ".txt" : ""));
        options.mleFileName = logFileNameField.getText();
    }

    public void setOptions(MarginalLikelihoodEstimationOptions options) {
        this.options = options;

        pathStepsField.setValue(options.pathSteps);
        chainLengthField.setValue(options.mleChainLength);
        logEveryField.setValue(options.mleLogEvery);

        logFileNameField.setText(options.mleFileName);

        treeWorkingPrior.setSelectedItem(options.choiceTreeWorkingPrior);
        if (options.choiceTreeWorkingPrior.equals("Product of exponential distributions") && options.performMLEGSS) {
            beautiOptions.logCoalescentEventsStatistic = true;
        } else {
            beautiOptions.logCoalescentEventsStatistic = false;
        }

        optionsPanel.validate();
        optionsPanel.repaint();
    }

    public void getOptions(MarginalLikelihoodEstimationOptions options) {
        options.pathSteps = pathStepsField.getValue();
        options.mleChainLength = chainLengthField.getValue();
        options.mleLogEvery = logEveryField.getValue();

        options.mleFileName = logFileNameField.getText();
        options.choiceTreeWorkingPrior = treeWorkingPrior.getSelectedItem().toString();
        if (options.choiceTreeWorkingPrior.equals("Product of exponential distributions") && options.performMLEGSS) {
            beautiOptions.logCoalescentEventsStatistic = true;
        } else {
            beautiOptions.logCoalescentEventsStatistic = false;
        }
    }

}
