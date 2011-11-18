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

DROP VIEW ComparisonSummaryView ;
CREATE VIEW ComparisonSummaryView AS
SELECT      concordant_callsets, COUNT(*)
FROM        SNPLevelView
GROUP BY    concordant_callsets
ORDER BY COUNT(*) DESC
;



SELECT      callset1, callset2, COUNT(*)
FROM        Comparison
WHERE       variable_id == 3 AND method_id ==1 AND value > 0.01
GROUP BY    callset1, callset2
LIMIT 10
;



DROP View SNPLevelView2 ;
CREATE VIEW IF NOT EXISTS SNPLevelView2 AS
SELECT V.rsid, V.chromosome, V.position, V.alleleA, V.alleleB, E.name AS comparison_method,
    ( SELECT COUNT( C.variant_id ) FROM Comparison C WHERE C.variant_id == V.id AND C.method_id = 1 AND C.variable_id == 3 AND C.value > 0.01 ) AS concordance_count
FROM    Variant V
INNER JOIN Entity E ON E.id == 2
GROUP BY V.chromosome, V.position, V.rsid
;
