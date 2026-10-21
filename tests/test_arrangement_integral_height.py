import json

from scripts.certify_arrangement_integral_height import main


def test_arrangement_integral_height_certificate(capsys):
    main()
    certificate = json.loads(capsys.readouterr().out)

    assert certificate["lambert"]["contracting_sup"] == "225/256"
    assert certificate["lambert"]["degree_bound"] == 92
    assert certificate["lehmer"]["trace_roots_in_minus2_2"] == 4
    assert certificate["lehmer"]["trace_roots_in_2_3"] == 1
    assert certificate["farey"]["log_volume"] == "n**2/4"
    assert certificate["farey"]["new_bad_prime_check"][-1] == 29
