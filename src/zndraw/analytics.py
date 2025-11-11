"""Prometheus metrics for ZnDraw analytics.

This module provides Prometheus metrics for monitoring ZnDraw application health
and usage.
"""

from prometheus_client import Gauge

# User metrics
connected_users = Gauge('connected_users', 'Number of currently connected users')

# Extension metrics - track registered/running extensions by category
running_modifiers = Gauge(
    'zndraw_running_modifiers',
    'Number of running modifier extensions',
    ['scope']  # Labels: 'room' or 'global'
)

running_selections = Gauge(
    'zndraw_running_selections',
    'Number of running selection extensions',
    ['scope']  # Labels: 'room' or 'global'
)

running_analysis = Gauge(
    'zndraw_running_analysis',
    'Number of running analysis extensions',
    ['scope']  # Labels: 'room' or 'global'
)
