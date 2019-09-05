**Updates in MathIOmica 1.2.1**
1. Source Code
	* Functions updates:
		* GOAnalysis and KEGGAnalysis default hypergeometric function updated to numeric evaluation to fix error in Mathematica 12.0.0 
2. Documentation
	* Updated (typographical errors and output updates):
		* GOAnalysis
		* KEGGAnalysis
		
**Updates in MathIOmica 1.2.0**
1. New Functions
	* ApplyBoxCoxTransformExtended
	* BoxCoxTransformExtended
2. Source Code
	* New private functions:
		* MaximumLikelihoodBoxCoxTransformBaseExtended
		* NormalLogLikelihoodBoxCoxBaseExtended
	* Functions updates:
		* ApplyBoxCoxTransformation now uses quiet mode, to fix error in Mathematica 12.0.0 from additional messaging.
		* GOAnalysisAssigner can address new GO Consortium resource GO file  
	* New File Updates:
		* MathIOmicaData:
			* goBasicObo.txt
			* humanIdentifierAssoc
			* humanGeneOntAssoc
			* humanGeneUCSCTable
3. Documentation
	* New:
		* ApplyBoxCoxTransformExtended
		* BoxCoxTransformExtended
	* Updated (typographical errors and output updates):
		* Autocorrelation
		* ApplyBoxCoxTransform
		* BoxCoxTransform
		* DataImporter
		* DataImporterDirect
		* FilteringFunction
		* FilterMissing
		* GetGeneDictionary
		* GOAnalysis
		* GOAnalysisAssigner
		* InverseAutocovariance
		* KEGGAnalysis
		* KEGGDictionary
		* KEGGPathwayVisuals
		* LombScargle
		* LowValueTag
		* MathIOmicaDynamicTranscriptome
		* MathIOmicaGuide
		* MathIOmicaOverview
		* MathIOmicaTutorial
		* MatrixClusters
		* MeasurementApplier
		* MSViewer
		* OBOGODictionary
		* OmicsObject
		* OmicsObjectCreator
		* QuantileEstimator
		* QuantileNormalization
		* Returner
		* SeriesApplier
		* StandardizeExtended
		* TimeSeriesClassification
		* TimeSeriesClusters
		* TimeSeriesSingleDendrogramHeatmap

**Updates in MathIOmica 1.1.3**
1. New Functions
	* None
2. Source Code
	* KEGG Pathway modification from http to https URLs.
3. Documentation
    * KEGGPathwayVisual Updated
    * KEGGAnalysisAssigner
    * KEGGDictionary
    * MathIOmica Overview
    * MathIOmica Dynamic Transcriptome
	
**Updates in MathIOmica 1.1.2**
1. New Functions
	* None
2. Source Code
	* Addressed performance issues when importing files, relating to Mathematica 11.2 updates of Import and incompatibility with specifying number of "HeaderLines". Internal code for importers updated to use Query for faster imports.
3. Documentation
	* None


**Updates in MathIOmica 1.1.1**
1. New Functions
	* None
2. Source Code
	* Addressed Issues when importing files relating to Mathematica 11.2 updates of Import and incompatibility with specifying number of "HeaderLines".
3. Documentation
	* MathIOmica Documentation was rebuilt on Mathematica Version 11.2.
	* When loading the Wolfram Help for the first time the MathIOmica help index in 1.1.0 is rebuilt. This is simply a warning with no errors, but an updated index version has been added.  

**Updates in MathIOmica 1.1.0**

1. New Functions 
	* OmicsObjectMerge
    * OmicsObjectPairMerge
2. Source Code
	* Updated code to reflect changes in behavior of DropKey in Mathematica Version 11.1.1
3. Documentation
	* Updated FilterMissing Documentation for typos
	* Updated TimeSeriesClassification for typos
	* Updated Dynamic Transcriptome Guide (Quantile Normalization specifications)
	* Updated MathIOmica Guide