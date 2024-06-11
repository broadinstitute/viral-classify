import unittest
from taxon_filter import calculate_chunking

class TestCalculateChunking(unittest.TestCase):
    def setUp(self):
        self.db_memory_estimate = 45  # GB

    def test_memory_exceeded(self):
        with self.assertRaises(ValueError):
            calculate_chunking(1000000, 16, 44, self.db_memory_estimate, 0.000002)  # Less than DB size

    def test_varying_read_sizes(self):
        # Test with small, medium, and large read sizes
        for read_size in [0.000001, 0.000002, 0.000005]:  # Different memory footprints
            chunk_sizes, threads_per_chunk = calculate_chunking(
                1000000, 16, 64, self.db_memory_estimate, read_size)
            self.assertTrue(all(size >= 20000 for size in chunk_sizes))  # Ensure all chunks meet minimum size

    def test_prevent_memory_oversaturation(self):
        # Simulate a scenario where each chunk could potentially use too much memory
        chunk_sizes, threads_per_chunk = calculate_chunking(
            1000000, 8, 64, 40, 0.000024)  # Set memory per read high
        total_memory_used = sum((cs * 0.000024 for cs in chunk_sizes)) + 40
        self.assertTrue(total_memory_used <= 64)  # Ensure total memory used does not exceed available memory

if __name__ == '__main__':
    unittest.main()
