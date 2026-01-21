import QtQuick
import JASP.Module

Description
{
	title		: qsTr("BSTS")
	description	: qsTr("This module offers a Bayesian take on linear Gaussian state space models suitable for   time series analysis.")
	icon     	: "bsts.png"
	hasWrappers	: true
	preloadData : false

	Analysis
	{
	    title: "Bayesian State Space Models"
	    func: "bayesianStateSpace"
		qml: 'bayesianStateSpace.qml'
	}
}
