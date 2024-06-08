import unittest
from taxon_filter import calculate_chunking


class TestCalculateChunking(unittest.TestCase):
    def setUp(self):
        self.db_memory_estimate = 45  # GB
        self.read_memory_per_read = 0.000002  # GB

    def test_basic_functional(self):
        chunk_sizes, threads_per_chunk = calculate_chunking(
            1000000, 16, 64, self.db_memory_estimate, self.read_memory_per_read)
        # Check if the total memory used is within limits
        # Assume each chunk uses its portion of the database memory plus reads memory
        total_memory_used = sum((cs * self.read_memory_per_read for cs in chunk_sizes)) + self.db_memory_estimate
        self.assertTrue(total_memory_used <= 64)
        self.assertGreaterEqual(threads_per_chunk, 1)

    def test_memory_constraint(self):
        chunk_sizes, threads_per_chunk = calculate_chunking(
            1000000, 16, 46, self.db_memory_estimate, self.read_memory_per_read)
        self.assertEqual(len(chunk_sizes), max(1, 1000000 // (int((46 - 45) / self.read_memory_per_read))))

    def test_high_read_count(self):
        chunk_sizes, threads_per_chunk = calculate_chunking(
            10000000, 32, 128, self.db_memory_estimate, self.read_memory_per_read)
        self.assertTrue(all(cs > 0 for cs in chunk_sizes))

    def test_minimal_thread_availability(self):
        chunk_sizes, threads_per_chunk = calculate_chunking(
            500000, 4, 64, self.db_memory_estimate, self.read_memory_per_read)
        self.assertEqual(threads_per_chunk, 1)

    def test_extremely_limited_memory(self):
        chunk_sizes, threads_per_chunk = calculate_chunking(
            500000, 16, 47, self.db_memory_estimate, self.read_memory_per_read)
        self.assertTrue(sum(chunk_sizes) >= 500000)  # Ensure all reads are processed

if __name__ == '__main__':
    unittest.main()
