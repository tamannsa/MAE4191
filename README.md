**MAE4191**

Master thesis code

The analyses in R produced 2,624 lines of code and can be accessed through GitHub. The OBD analyses were conducted in Excel using equation (3) found in the main master thesis document. The main workflow is outlined below.

**1** Packages 

**2** Load data 

**3** Variable selection 

**4** Missingness diagnostics 

	**4.1** Investigating missingness in scales 
 
	**4.2** Investigating missingness in items
 
	**4.3** Visualization and statistical tests 
 

**5** Multiple imputation
	**5.1** Identify and select auxiliary variables 
	**5.2** Split dataset with one plausible value each 
	**5.3** Split dataset by country 
	**5.4** Run test imputation 
		**5.4.1** Test imputation diagnostics
	**5.5** Run full imputation
		**5.5.1** Full imputation restructuring 
		**5.5.2** Full imputation diagnostics

**6** Multicollinearity diagnostics 

**7** Main analysis
	**7.1** Create survey design 
	**7.2** Regression analysis
	**7.3** Model comparison
	**7.4** Calculate R2
	**7.5** Linear model assumptions
 		**7.5.1** Normality assumption 
   		**7.5.2** Homoscedasticity assumption 
     		**7.5.3** Linearity assumption 

**8** Investigating outliers
	**8.1** Before imputation
	**8.2** After imputation

**9** Statistics 
	**9.1** Mean and standard deviation
	**9.2** Skewness and kurtosis 
	**9.3** Correlation
	**9.4** Sample sizes
