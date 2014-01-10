

CREATE TABLE SimilarSequences_prots (
 QUERY_ID                 VARCHAR(60),
 SUBJECT_ID               VARCHAR(60),
 QUERY_TAXON_ID           VARCHAR(40),
 SUBJECT_TAXON_ID         VARCHAR(40),
 EVALUE_MANT              FLOAT,
 EVALUE_EXP               INT,
 PERCENT_IDENTITY         FLOAT,
 PERCENT_MATCH            FLOAT
) 


CREATE INDEX ss_qtaxexp_ix_prots
ON SimilarSequences_prots(query_id, subject_taxon_id,
evalue_exp, evalue_mant,
query_taxon_id, subject_id)  


CREATE INDEX ss_seqs_ix_prots
ON SimilarSequences_prots(query_id, subject_id,
evalue_exp, evalue_mant, percent_match)  


CREATE TABLE InParalog_prots (
 SEQUENCE_ID_A           VARCHAR(60),
 SEQUENCE_ID_B           VARCHAR(60),
 TAXON_ID                VARCHAR(40),
 UNNORMALIZED_SCORE      FLOAT,
 NORMALIZED_SCORE        FLOAT
)


CREATE INDEX inparalog_seqa_ix_prots
ON InParalog_prots(sequence_id_a)  


CREATE INDEX inparalog_seqb_ix_prots
ON InParalog_prots(sequence_id_b)  

 
CREATE TABLE Ortholog_prots (
 SEQUENCE_ID_A           VARCHAR(60),
 SEQUENCE_ID_B           VARCHAR(60),
 TAXON_ID_A              VARCHAR(40),
 TAXON_ID_B              VARCHAR(40),
 UNNORMALIZED_SCORE      FLOAT,
 NORMALIZED_SCORE        FLOAT
)


CREATE INDEX ortholog_seq_a_ix_prots
ON Ortholog_prots(sequence_id_a)  


CREATE INDEX ortholog_seq_b_ix_prots
ON Ortholog_prots (sequence_id_b)  


CREATE TABLE CoOrtholog_prots (
 SEQUENCE_ID_A           VARCHAR(60),
 SEQUENCE_ID_B           VARCHAR(60),
 TAXON_ID_A              VARCHAR(40),
 TAXON_ID_B              VARCHAR(40),
 UNNORMALIZED_SCORE      FLOAT,
 NORMALIZED_SCORE        FLOAT
)


CREATE INDEX coortholog_seq_a_ix_prots
ON CoOrtholog_prots(sequence_id_a)  


CREATE INDEX coortholog_seq_b_ix_prots
ON CoOrtholog_prots (sequence_id_b)  


CREATE OR REPLACE VIEW InterTaxonMatch_prots
	AS SELECT ss.query_id, ss.subject_id, ss.subject_taxon_id,
	ss.evalue_mant, ss.evalue_exp
	FROM SimilarSequences_prots ss
	WHERE ss.subject_taxon_id != ss.query_taxon_id
