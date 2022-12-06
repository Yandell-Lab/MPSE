with nicus as (
	SELECT xwalk.dw_pid, xwalk.pat_id, pat.birth_date
	FROM visit_dm.visit_detail visit
	JOIN etl_process.dwid_crosswalk xwalk
	ON visit_dwid = dw_vid
	JOIN visit_dm.patient pat
	ON pat.pat_id = xwalk.pat_id
	AND visit.nicu_unit_flg_dwid = 675863
	WHERE pat.pat_id = '{{ placeholder }}'
)

SELECT nicus.pat_id,
nicus.birth_date, clobs.rpt_id,
clobs.proc_dt, clobs.rpt_type_cd, clobs.text
FROM clindata.general_text_clobs clobs
JOIN nicus
ON nicus.pat_id = clobs.pat_id
AND MONTHS_BETWEEN(clobs.proc_dt, nicus.birth_date) <= 1
