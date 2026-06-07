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
