import json
import pytest
import g2ccx_paths as gp


@pytest.fixture(autouse=True)
def clean(monkeypatch, tmp_path):
    monkeypatch.delenv(gp.ENV_SOLVER, raising=False)
    monkeypatch.chdir(tmp_path)  # no local file here
    monkeypatch.setattr(gp, "_local_file_solver", lambda: (None, None))
    monkeypatch.setattr(gp.shutil, "which", lambda name: None)


def exe(tmp_path, name="ccx.exe"):
    p = tmp_path / name
    p.write_text("")
    return str(p)


def test_solver_priority_config_env_local_path(tmp_path, monkeypatch):
    a, b, c, d = (exe(tmp_path, n) for n in ("a", "b", "c", "d"))
    monkeypatch.setattr(gp.shutil, "which", lambda name: d if name == "ccx" else None)
    assert gp.resolve_solver() == d
    monkeypatch.setattr(gp, "_local_file_solver", lambda: (c, None))
    assert gp.resolve_solver() == c
    monkeypatch.setenv(gp.ENV_SOLVER, b)
    assert gp.resolve_solver() == b
    assert gp.resolve_solver(a) == a


def test_missing_solver_explains_the_options(tmp_path, monkeypatch):
    monkeypatch.setenv(gp.ENV_SOLVER, str(tmp_path / "nope"))
    with pytest.raises(FileNotFoundError, match="G2CCX_SOLVER.*not a file|not a file"):
        gp.resolve_solver()


def test_output_timestamp_create_and_find(tmp_path):
    pattern = str(tmp_path / "{timestamp}_run")
    assert gp.resolve_output(str(tmp_path / "plain"), True) == str(tmp_path / "plain")
    new = gp.resolve_output(pattern, True)
    assert new.endswith("_run") and len(new.rsplit("\\", 1)[-1].rsplit("/", 1)[-1]) == 12 + len("_run")
    with pytest.raises(FileNotFoundError):
        gp.resolve_output(pattern, False)
    (tmp_path / "202601010000_run").mkdir()
    (tmp_path / "202602020000_run").mkdir()
    (tmp_path / "202603030000_other").mkdir()
    assert gp.resolve_output(pattern, False).endswith("202602020000_run")
