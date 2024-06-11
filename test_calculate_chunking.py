import unittest
from taxon_filter import calculate_chunking

class TestCalculateChunking(unittest.TestCase):
    def setUp(self):
        self.db_memory_estimate = 45  # GB

    def test_memory_constrained_environment(self):
        # Settings where the database size is near the total available memory
        number_of_reads = 1000000
        total_threads = 16
        max_memory = 46  # GB, slightly more than the database estimate
        db_memory_estimate = 45  # GB, almost the total available memory
        read_memory_per_read = 0.000002  # GB per read

        chunk_sizes, threads_per_chunk = calculate_chunking(
            number_of_reads, total_threads, max_memory, db_memory_estimate, read_memory_per_read
        )

        # Expect minimal number of chunks and all threads to be used for likely a single chunk
        expected_chunks = 1
        expected_threads = total_threads

        # Assertions to check if the function behaves as expected
        self.assertEqual(len(chunk_sizes), expected_chunks, "Function should minimize number of chunks under tight memory conditions")
        self.assertEqual(threads_per_chunk, expected_threads, "Function should use all available threads for the minimal number of chunks")

        # Additionally, ensure that the memory total is within the bounds we gave it 
        total_memory_used = sum(size * read_memory_per_read for size in chunk_sizes) + db_memory_estimate
        self.assertTrue(total_memory_used <= max_memory, "Total memory usage should not exceed available memory")

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
