"""注釈の出所と候補選択を保つ回帰テスト。"""
import copy

import pytest

from metabolomix.core import path_resolvers


def test_catalog_annotation_selects_unknown_arf_without_mutating_input():
    spots = [{'MasterAlignmentID': 1, 'Name': 'Unknown', 'HeightAverage': 100},
             {'MasterAlignmentID': 2, 'Name': 'PC 34:1', 'HeightAverage': 200}]
    original = copy.deepcopy(spots)
    catalog = [{'MasterAlignmentID': 1, 'Name': 'SL 33:0;O', 'Ontology': 'SL'},
               {'MasterAlignmentID': 2, 'Name': 'LPC 18:0', 'Ontology': 'LPC'}]
    selected = path_resolvers._filter_arf_spots(spots, annotation_keyword='SL', catalog=catalog)
    assert [p['MasterAlignmentID'] for p in selected] == [1]
    assert selected[0]['Name'] == 'SL 33:0;O'
    assert selected[0]['annotation_source'] == 'arf2'
    assert selected[0]['arf_name'] == 'Unknown'
    assert selected[0]['annotation_conflict'] is True
    assert spots == original


def test_ontology_exact_match_and_id_intersection_do_not_include_lpc():
    spots = [{'MasterAlignmentID': i, 'Name': name} for i, name in enumerate(['PC 34:1', 'LPC 18:0', 'PC 36:2'])]
    catalog = [{'MasterAlignmentID': i, 'Name': name, 'Ontology': ontology}
               for i, name, ontology in [(0,'PC 34:1','PC'),(1,'LPC 18:0','LPC'),(2,'PC 36:2','PC')]]
    selected = path_resolvers._filter_arf_spots(spots, catalog=catalog, ontologies=['PC'], spot_ids=[0,1])
    assert [p['MasterAlignmentID'] for p in selected] == [0]


@pytest.mark.parametrize('kwargs', [{'spot_ids':[]}, {'spot_ids':[True]}, {'spot_ids':[-1]},
                                  {'spot_ids':[1.5]}, {'ontologies':[]}, {'ontologies':['']}, {'ontologies':[None]}])
def test_invalid_selectors_do_not_silently_broaden_selection(kwargs):
    with pytest.raises(ValueError):
        path_resolvers._filter_arf_spots([{'MasterAlignmentID':1}], **kwargs)


def test_catalog_duplicate_id_is_not_silently_overwritten():
    with pytest.raises(ValueError):
        path_resolvers._filter_arf_spots([{'MasterAlignmentID':1}], catalog=[
            {'MasterAlignmentID':1,'Name':'PC'}, {'MasterAlignmentID':1,'Name':'SL'}])


def test_unknown_catalog_name_falls_back_to_arf_and_no_match_stays_empty():
    spots = [{'MasterAlignmentID':1,'Name':'SL 33:0;O'}]
    assert path_resolvers._filter_arf_spots(spots, annotation_keyword='SL',
        catalog=[{'MasterAlignmentID':1,'Name':'Unknown'}])[0]['annotation_source'] == 'arf'
    assert path_resolvers._filter_arf_spots(spots, spot_ids=[999]) == []


@pytest.mark.parametrize('value', [True, 1.0, '1'])
def test_mcp_rejects_coerced_spot_ids(value):
    import server
    from pydantic import ValidationError
    model = server.mcp._tool_manager.get_tool('arf_parser').fn_metadata.arg_model
    with pytest.raises(ValidationError):
        model.model_validate({'spot_ids':[value]})
