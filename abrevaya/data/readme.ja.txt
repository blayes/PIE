J. Abrevaya, "Estimating the Effect of Smoking on Birth Outcomes Using a
Matched Panel Data Approach," Journal of Applied Econometrics, Vol. 21,
No. 4, 2006, pp. 489-519.

There is one data file, which is an ASCII file in DOS format. This file,
birpanel.txt, is zipped in the file birpanel.zip. Unix users should use
"unzip -a".

These data represent a subsample of the 1990-1998 Natality Data Sets
which are publicly available through the National Center for Health
Statistics. The subsample was constructed through a matching algorithm
that is described in detail in Section 2 of the paper.

The file birpanel.txt contains the data corresponding to "Matched Panel
#3" in the paper.  The file contains 296,218 observations, with 141,929
distinct mothers (identified by the variable "momid3" described below). 
The variables (whose names appear in the first line of the text file)
are as follows:

momid3             identification number of the mother
idx                index number of a mother's birth (starts with 1 and goes
                   to at least 2 for each mother); momid3 and idx together
                   describe the panel structure of the data
stateres           code (1 through 51) for state of residence; see the NCHS
                   (National Center for Health Statistics) website for codes
dmage              age of mother (in years)
dmeduc             education of mother (in years)
mplbir             state code (1 through 51) for mother's state of birth
nlbnl              number of live births now living
gestat             length of gestation (in weeks)
dbirwt             birthweight (in grams)
cigar              number of cigarettes smoked per day (99=unknown)
smoke              indicator variable for smoking status (1=smoker,
		   0=nonsmoker)
male               indicator variable for baby gender (1=male, 0=female)
year               year of birth (0=1990, ..., 8=1998)
married            indicator variable for marital status (1=married,
		   0=unmarried)
hsgrad             high-school graduate indicator (constructed from dmeduc)
somecoll           some-college indicator (constructed from dmeduc)
collgrad           college-graduate indicator (constructed from dmeduc)
agesq              age of mother squared
black              indicator variable for black race (1=black, 0=white)
adeqcode2          indicator that Kessner index = 2
adeqcode3          indicator that Kessner index = 3
novisit            indicator that no prenatal visit occurred
pretri2            indicator that first prenatal visit occurred in 2nd
		   trimester
pretri3            indicator that first prenatal visit occurred in 3nd
		   trimester
proxy_exists       indicator variable for whether or not the proxy variable
		   for a correct match exists for this observation: if equal
		   to 1, the following variable contains the value of the
		   proxy; if equal to 0, the following variable contains
		   "proxyhat" (see below)
proxy_or_proxyhat  if proxy_exists = 1, this contains a proxy variable for
		   correct match (1 if interval since last live birth
		   matched the actual birthdate of the previous birth); if
		   proxy_exists = 0, this contains the predicted probability
		   that proxy variable is equal to one (see Section 4.2 of
		   the paper; this variable was constructed from a probit
		   regression using a mother's first and second births; the
		   probit specification included birthorder dummies,
		   education-level dummies, age, age squared, race, and
		   marital status)

The results in the paper using the last two variables (proxy or predicted
proxy) are based only on a mother's first two births (idx<=2).  The other
results (that do not utilize the proxy information) are based on all births.
