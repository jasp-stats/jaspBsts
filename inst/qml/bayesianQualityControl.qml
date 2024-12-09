import QtQuick
import QtQuick.Layouts
import JASP
import JASP.Controls

Form
{
	columns: 2

	VariablesForm
	{

		AvailableVariablesList	{ name: "allVariablesList" }
		AssignedVariablesList	{ name: "measurement";		title: qsTr("Measurement");		allowedColumns: 	["scale"];	singleVariable: true}
		AssignedVariablesList	{ name: "time";				title: qsTr("Time");			allowedColumns: 	["scale"];	singleVariable: true}

		Group
		{
			title: qsTr("Control Limits")
			columns: 1


			DoubleField
			{
				name: 			"controlLimitsLower"
				id:				controlLimitsLower
				label:			qsTr("Lower")
				defaultValue:	-1
				min:			-Infinity
				max:			controlLimitsUpper.value
			}

			DoubleField
			{
				name: 			"controlLimitsUpper"
				id:				controlLimitsUpper
				label:			qsTr("Upper")
				defaultValue:	1
				min:			controlLimitsLower.value
				max:			Infinity
			}
		}
	}

	Group
	{
		title:	qsTr("State Plots")

		CheckBox
		{
			name:		"statePlotsMean"
			label:		qsTr("Mean")
			checked:	true
		}

		CheckBox
		{
			name:		"statePlotsStandardDeviation"
			label:		qsTr("Standard deviation")
			checked:	true
		}

		CheckBox
		{
			name:		"statePlotsCpk"
			label:		qsTr("Cₚₖ")
			checked:	true
		}
	}

	Group
	{
		title:	qsTr("Posterior Plots")

		CheckBox
		{
			name:	"posteriorPlotsMean"
			label:	qsTr("Mean")	
		}

		CheckBox
		{
			name:	"posteriorPlotsStandardDeviation"
			label:	qsTr("Standard deviation")	
		}

		CheckBox
		{
			name:	"posteriorPlotsCpk"
			label:	qsTr("Cₚₖ")	
		}

		IntegerField
		{
			name:			"posteriorPlotsAtState"
			label:			qsTr("At state")
			defaultValue:	1
		}
	}

	Section
	{
		title:	qsTr("Advanced")

		Group
		{
			title: 		qsTr("MCMC")

			IntegerField
			{
				name:			"advancedMcmcAdaptation"
				label:			qsTr("Adaptation")
				defaultValue:	500
				min:			100
				fieldWidth:		55 * preferencesModel.uiScale
			}
			IntegerField
			{
				name:			"advancedMcmcBurnin"
				label:			qsTr("Burnin")
				defaultValue:	2000
				min:			100
				fieldWidth:		55 * preferencesModel.uiScale
			}
			IntegerField
			{
				name:			"advancedMcmcSamples"
				label:			qsTr("Samples")
				defaultValue:	5000
				min:			100
				fieldWidth:		55 * preferencesModel.uiScale
			}
			IntegerField
			{
				name:			"advancedMcmcChains"
				label:			qsTr("Chains")
				defaultValue:	3
				min:			1
				fieldWidth:		55 * preferencesModel.uiScale
			}
			IntegerField
			{
				name:			"advancedMcmcThin"
				label:			qsTr("Thin")
				defaultValue:	1
				min:			1
				fieldWidth:		55 * preferencesModel.uiScale
			}
		}

		IntegerField
		{
			name:			"advancedStateAggregation"
			label:			qsTr("State aggregation")
			defaultValue:	20
		}

	}
}
