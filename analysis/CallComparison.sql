/* Convert to chromosome/position/rsid index */
DROP INDEX Variant_index ;
CREATE INDEX IF NOT EXISTS Variant_index ON Variant( chromosome, position, rsid ) ;
DROP INDEX ComparisonIndex ;
CREATE INDEX IF NOT EXISTS ComparisonIndex ON Comparison( variant_id, method_id, variable_id )

/* Create views of the data */
CREATE VIEW IF NOT EXISTS ComparisonOverview AS
SELECT V.rsid, V.chromosome, V.position, E.name AS comparison_method, MIN( C.value ) AS min_pvalue, MAX( C.value ) AS max_pvalue
FROM    Variant V
INNER JOIN Comparison C
    ON  C.variant_id == V.id AND C.method_id == 1 AND C.variable_id == 3
INNER JOIN Entity E
    ON E.id == C.method_id
GROUP BY V.chromosome, V.position, V.rsid
;

DROP VIEW ComparisonOverview ;
CREATE VIEW IF NOT EXISTS ComparisonOverview AS
SELECT V.rsid, V.chromosome, V.position, E.name AS comparison_method, C.callset1, C.callset2, C.value AS pvalue
FROM    Variant V
INNER JOIN Comparison C
    ON  C.variant_id == V.id AND C.method_id == 1 AND C.variable_id == 3
INNER JOIN Entity E
    ON E.id == C.method_id
;

DROP View SNPLevelView ;
CREATE VIEW IF NOT EXISTS SNPLevelView AS
SELECT V.rsid, V.chromosome, V.position, V.alleleA, V.alleleB, COUNT( C.variant_id ) AS concordance_count, GROUP_CONCAT( C.callset1 || ":" || C.callset2, "," ) AS concordant_callsets
FROM    Variant V
LEFT OUTER JOIN Comparison C
    ON  C.variant_id == V.id AND C.method_id == 1 AND C.variable_id == 3 AND C.value > 0.001
GROUP BY V.chromosome, V.position, V.rsid
;

CREATE VIEW IF NOT EXISTS SNPLevelView_1E04 AS
SELECT V.rsid, V.chromosome, V.position, V.alleleA, V.alleleB, COUNT( C.variant_id ) AS concordance_count, GROUP_CONCAT( C.callset1 || ":" || C.callset2, "," ) AS concordant_callsets
FROM    Variant V
LEFT OUTER JOIN Comparison C
    ON  C.variant_id == V.id AND C.method_id == 1 AND C.variable_id == 3 AND C.value > 0.0001
GROUP BY V.chromosome, V.position, V.rsid
;

CREATE VIEW IF NOT EXISTS SNPLevelView_1E05 AS
SELECT V.rsid, V.chromosome, V.position, V.alleleA, V.alleleB, COUNT( C.variant_id ) AS concordance_count, GROUP_CONCAT( C.callset1 || ":" || C.callset2, "," ) AS concordant_callsets
FROM    Variant V
LEFT OUTER JOIN Comparison C
    ON  C.variant_id == V.id AND C.method_id == 1 AND C.variable_id == 3 AND C.value > 0.00001
GROUP BY V.chromosome, V.position, V.rsid
;

CREATE VIEW IF NOT EXISTS SNPLevelView_1E06 AS
SELECT V.rsid, V.chromosome, V.position, V.alleleA, V.alleleB, COUNT( C.variant_id ) AS concordance_count, GROUP_CONCAT( C.callset1 || ":" || C.callset2, "," ) AS concordant_callsets
FROM    Variant V
LEFT OUTER JOIN Comparison C
    ON  C.variant_id == V.id AND C.method_id == 1 AND C.variable_id == 3 AND C.value > 0.000001
GROUP BY V.chromosome, V.position, V.rsid
;

DROP VIEW ComparisonSummaryView ;
CREATE VIEW ComparisonSummaryView AS
SELECT      concordant_callsets, concordance_count, COUNT(*)
FROM        SNPLevelView
GROUP BY    concordant_callsets
ORDER BY COUNT(*) DESC
;

CREATE VIEW ComparisonSummaryView1E04 AS
SELECT      concordant_callsets, concordance_count, COUNT(*)
FROM        SNPLevelView_1E04
GROUP BY    concordant_callsets
ORDER BY COUNT(*) DESC
;

CREATE VIEW ComparisonSummaryView1E05 AS
SELECT      concordant_callsets, concordance_count, COUNT(*)
FROM        SNPLevelView_1E05
GROUP BY    concordant_callsets
ORDER BY COUNT(*) DESC
;

CREATE VIEW ComparisonSummaryView1E06 AS
SELECT      concordant_callsets, concordance_count, COUNT(*)
FROM        SNPLevelView_1E06
GROUP BY    concordant_callsets
ORDER BY COUNT(*) DESC
;
