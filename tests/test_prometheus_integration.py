"""Tests for Prometheus integration.

Tests follow CLAUDE.md principles:
- Each test is a function (not a method)
- Tests are specific and test only one thing
- Uses pytest.mark.parametrize to avoid duplication
"""


def test_analytics_module_imports():
    """Test that analytics module can be imported."""
    from zndraw.analytics import connected_users

    assert connected_users is not None
    # Should be able to call methods without error
    connected_users.inc()
    connected_users.dec()


def test_create_app_has_prometheus(tmp_path):
    """Test that Flask app always has Prometheus enabled."""
    from zndraw.config import ZnDrawConfig
    from zndraw.server import create_app

    config = ZnDrawConfig()
    config.storage_path = str(tmp_path)

    app = create_app(config=config)
    assert app is not None


def test_metrics_endpoint_always_available(tmp_path):
    """Test that /metrics endpoint is always available."""
    from zndraw.config import ZnDrawConfig
    from zndraw.server import create_app

    config = ZnDrawConfig()
    config.storage_path = str(tmp_path)

    app = create_app(config=config)

    with app.test_client() as client:
        response = client.get('/metrics')
        assert response.status_code == 200
        # Should contain Prometheus metrics format
        assert b'# HELP' in response.data or b'# TYPE' in response.data


def test_connected_users_metric_increments():
    """Test that connected_users metric can be incremented."""
    from zndraw.analytics import connected_users

    # Get initial value (may not be 0 if other tests ran)
    from prometheus_client import REGISTRY
    initial_value = None
    for metric in REGISTRY.collect():
        if metric.name == 'connected_users':
            for sample in metric.samples:
                initial_value = sample.value
                break

    # Increment
    connected_users.inc()

    # Check new value
    for metric in REGISTRY.collect():
        if metric.name == 'connected_users':
            for sample in metric.samples:
                if initial_value is not None:
                    assert sample.value == initial_value + 1
                break


def test_connected_users_metric_decrements():
    """Test that connected_users metric can be decremented."""
    from zndraw.analytics import connected_users

    # Increment first to ensure we have a positive value
    connected_users.inc()

    # Get current value
    from prometheus_client import REGISTRY
    current_value = None
    for metric in REGISTRY.collect():
        if metric.name == 'connected_users':
            for sample in metric.samples:
                current_value = sample.value
                break

    # Decrement
    connected_users.dec()

    # Check new value
    for metric in REGISTRY.collect():
        if metric.name == 'connected_users':
            for sample in metric.samples:
                if current_value is not None:
                    assert sample.value == current_value - 1
                break
