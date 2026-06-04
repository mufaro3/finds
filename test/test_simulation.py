def test_serialize_to_file():
    r"""
    For an input shape of :math:`(N,6)`:
    
    1. Each write should create only one row.
    2. Each row should have exactly :math:`6N+1` columns.
    3. The deserialized version of the data should be identical
       to the data before serialization. (We should be able to
       fully retrieve the data afterward.)
    """
    pass

def test_calculate_update_rk4():
    r"""
    For an input shape of :math:`(N,6)`:
    
    1. The output shape should always be the same as the input shape.
    2. For a time step of 0, the output should be identical to the input.
    3. Smaller time steps should produce less change than bigger time steps,
       or :math:`|| \mathbf{X} - \mathbf{X}_{\delta t_1} || \le || \mathbf{X}
       - \mathbf{X}_{\delta t_2} ||` where :math:`\delta t_1 < \delta t_2`
    4. The Eulerian Limit: for small :math:`\delta t`, Eulerian approximation
       and Runge-Kutta should be roughly equal.
    
       .. math::
           \mathbf{X}(t + \delta t)_{\text{RK4}} \approx \mathbf{X} +
           f(\mathbf{X}) \delta t

    5. The orientation vector should come out normalized.
    6. The calculation should be deterministic, so running the same
       function twice should produce identical results.
    """
    pass

def test_generate_animation_from_file():
    pass

def test_perform_simulation():
    pass
