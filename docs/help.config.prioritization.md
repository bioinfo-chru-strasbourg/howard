# HOWARD Prioritization configuration

## Prioritization Configuration Guide

Prioritization algorithm uses profiles to flag variants (as passed or filtered), calculate a prioritization score, and automatically generate a comment for each variants (example: 'polymorphism identified in dbSNP. associated to Lung Cancer. Found in ClinVar database'). Prioritization profiles are defined in a configuration file in JSON format. A profile is defined as a list of annotation/value, using wildcards and comparison options (contains, lower than, greater than, equal...). Annotations fields may be quality values (usually from callers, such as 'DP') or other annotations fields provided by annotations tools, such as HOWARD itself (example: COSMIC, Clinvar, 1000genomes, PolyPhen, SIFT). 

See [HOWARD Help Prioritization tool](help.md#prioritization-tool) tool for more information.

## Example

This example describes the prioritization profile 'default' that uses 2 fields 'DP' ('Read Depth') and 'CLNSIG' ('ClinVar Significance') with specific criteria.

For 'DP' field, 2 filters are applied: if 'DP' is greater than or equal to '50', score is '5' and flag is 'PASS', if 'DP' is lower than '50', score is '0' and flag is 'FILTERED'. It means that if 'Read depth' is lower than 50 for a variant, it will be filtered with a bad score. Otherwise, the variant will have a better score (but can be filtered because of other filters).

For field 'CLNSIG', 2 filters are applied: if it is equal to 'pathogenic', score is '15' and flag is 'PASS', and if it is equal to 'non-pathogenic', score is '-100' and flag is 'FILTERED'. Thus, the variant will be well scored if it is pathogenic, but filtered with a bad score if it is non pathogenic (for Clinvar annotation).

```json
{
    "default": {
        "DP": [
            {
                "type": "gte",
                "value": "50",
                "score": 5,
                "flag": "PASS",
                "comment": [
                    "DP higher than 50"
                ]
            },
            {
                "type": "lt",
                "value": "50",
                "score": 0,
                "flag": "FILTERED",
                "comment": [
                    "DP lower than 50"
                ]
            }
        ],
        "CLNSIG": [
            {
                "type": "equals",
                "value": "pathogenic",
                "score": 15,
                "flag": "PASS",
                "comment": [
                    "Described on CLINVAR database as pathogenic"
                ]
            },
            {
                "type": "equals",
                "value": "non-pathogenic",
                "score": -100,
                "flag": "FILTERED",
                "comment": [
                    "Described on CLINVAR database as non-pathogenic"
                ]
            }
        ]
    }
}

```

## Prioritization fields explanation

### Prioritization Profile

A prioritization profile contains filters applied to specific fields to prioritize or filter variants based on defined criteria. Each profile is identified by a unique name (e.g. 'default') , and defines filters applied to a specific field (e.g. 'DP', 'CLNSIG'). Each filter consists of criteria such as type of test (e.g. greater than, contains), threshold value or substring, and return a score, a flag, and a comment if the test is valid.

### Filter criteria

#### type
  - *Type:* String
  - *Description:* Specifies the type of test to apply to the field. It can be 
     - 'gt' for 'greater than' ('>')
     - 'gte' for 'greater than or equal to' ('>=')
     - 'lt' for 'lower than' ('<')
     - 'lte' for 'lower than or equal to' ('<=')
     - 'equals' for exact match ('='), for integer or string comparison
     - 'contains' for exact match

#### value
  - *Type:* String or Integer
  - *Description:* Specifies the threshold value or string/substring to match within the field, depending on the type of test specified.

### Filter result

#### score
  - *Type:* Integer
  - *Description:* Assigns a score to variants that pass the filter. The score is used for prioritization or ranking purposes. Depending on the prioritization calculation mode, this score will be incremented ('HOWARD' mode) or compared if it is the max ('VaRank' mode).

#### flag
  - *Type:* String
  - *Description:* Assigns a flag (either 'PASS' or 'FILTERED') to variants that pass the filter. The flag provides additional information or categorization about the variant. Flag 'FILTERED' is prior to 'PASS' to calculate the global prioritization profile flag. 

- **comment:**
  - *Type:* Array of Strings
  - *Description:* Provides a comment or explanation for the filter criteria. It helps users understand the rationale behind applying the filter.

