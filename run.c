#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX 0x100
#define METRIC sum // sum, minus, max
// #define BOUND
const char charset[] = "01";

// Based on the last character, count the number of times a sequence repeats itself, multiplied by the sequence length.
int count_sum(char *sequence, int base)
{
    int count = 0;
    for (int offset = base - 1; offset >= 0; offset--)
    {
        for (int i = offset, j = base; i >= 0 && sequence[i] == sequence[j]; i--, j--)
        {
            count++;
        }
    }
    return count;
}

// Based on the last character, count the number of times a sequence repeats itself, multiplied by the sequence length minus one.
int count_minus(char *sequence, int base)
{
    int sum = 0;
    for (int offset = base - 1; offset >= 0; offset--)
    {
        int count = 0;
        for (int i = offset, j = base; i >= 0 && sequence[i] == sequence[j]; i--, j--)
        {
            count++;
        }
        if (count > 1)
        {
            sum += count - 1;
        }
    }
    return sum;
}

// Based on the last character, return the length of the longest sequence repeats itself.
int count_max(char *sequence, int base)
{
    int max = 0;
    for (int offset = base - 1; offset >= 0; offset--)
    {
        int count = 0;
        for (int i = offset, j = base; i >= 0 && sequence[i] == sequence[j]; i--, j--)
        {
            count++;
        }
        if (max < count)
        {
            max = count;
        }
    }
    return max;
}

#define CONCAT_(a, b) a##b
#define CONCAT(a, b) CONCAT_(a, b)
#define count(...) CONCAT(count_, METRIC)(__VA_ARGS__)

// Get the total repeat count by summing the count of each item.
long long accumulate_trunc(char *sequence, int length)
{
    long long sum = 0;
    for (int i = 0; i < length; i++)
    {
        sum += count(sequence, i);
    }
    return sum;
}

// Parse a string into a sequence of character indices.
int parse(const char *string, char *sequence, int length)
{
    int base = length < (int)strlen(string) ? length : (int)strlen(string);
    for (int i = 0; i < base; i++)
    {
        sequence[i] = strchr(charset, string[i]) - charset;
    }
    return base;
}

// Serialize a sequence of character indices into a string.
int serialize(const char *sequence, int length, char *string)
{
    for (int i = 0; i < length; i++)
    {
        string[i] = charset[(int)sequence[i]];
    }
    string[length] = '\0';
    return length;
}

// Compare two sequences.
int compare(char *sa, char *sb, int length)
{
    for (int i = 0; i < length; i++)
    {
        if (sa[i] != sb[i])
        {
            return sa[i] - sb[i];
        }
    }
    return 0;
}

// Print a sequence to a file pointer, e.g. stdout.
void print(char *sequence, int length, FILE *file)
{
    for (int i = 0; i < length; i++)
    {
        fprintf(file, "%c", charset[(int)sequence[i]]);
    }
}

#define STR_(s) #s
#define STR(s) STR_(s)

// Generate a file name based on the metric and character set.
void name_ref(char *buffer, int length)
{
    snprintf(buffer, length, "reference/%s_%s.txt", STR(METRIC), charset);
}

// List of character indices that make up a sequence.
typedef struct
{
    int repeat;
    char *sequence;
} sequence_t;

// List of sequences that make up a reference collection to compare generated sequences agains.
typedef struct
{
    int length;
    float offset;
    float divide;
    float add;
    sequence_t *sequences;
} reference_t;

// Read a reference collection from a file.
void read_ref(reference_t *reference)
{
    char buffer[MAX + 1] = {0};
    name_ref(buffer, MAX);
    FILE *file = fopen(buffer, "r");
    if (!file)
    {
        return;
    }
    int count = 1;
    fscanf(file, "%d %f %f %f", &count, &(reference->offset), &(reference->divide), &(reference->add));
    sequence_t *sequences = malloc(count * sizeof(sequence_t));
    sequences[0].repeat = 0;
    sequences[0].sequence = "\0";
    for (int length = 2; length <= count; length++)
    {
        fscanf(file, "%d %s", &(sequences[length - 1].repeat), buffer);
        sequences[length - 1].sequence = malloc(length * sizeof(char));
        parse(buffer, sequences[length - 1].sequence, length);
    }
    fclose(file);
    reference->length = count;
    reference->sequences = sequences;
}

// Write a reference collection to a file.
void write_ref(reference_t reference)
{
    char buffer[MAX + 1] = {0};
    name_ref(buffer, MAX);
    FILE *file = fopen(buffer, "w");
    fprintf(file, "%d %.2f %.2f %.2f\n", reference.length, reference.offset, reference.divide, reference.add);
    sequence_t *sequences = reference.sequences;
    for (int length = 2; length <= reference.length; length++)
    {
        serialize(sequences[length - 1].sequence, length, buffer);
        fprintf(file, "%d %s\n", sequences[length - 1].repeat, buffer);
    }
    fclose(file);
}

// Merge a sequence into a reference collection and write if improved.
void merge_ref(char *sequence, int length, reference_t *reference)
{
    if (length <= reference->length)
    {
        sequence_t *sequences = reference->sequences;
        int repeat = accumulate_trunc(sequence, length);
        if (repeat < sequences[length - 1].repeat)
        {
            sequences[length - 1].repeat = repeat;
            memcpy(sequences[length - 1].sequence, sequence, length * sizeof(char));
            write_ref(*reference);
        }
        else if (repeat == sequences[length - 1].repeat)
        {
            int comp = compare(sequence, sequences[length - 1].sequence, length);
            if (comp < 0)
            {
                memcpy(sequences[length - 1].sequence, sequence, length * sizeof(char));
                write_ref(*reference);
            }
        }
    }
    else if (length == reference->length + 1)
    {
        sequence_t *sequences = realloc(reference->sequences, (length + 1) * sizeof(sequence_t));
        sequences[length - 1].repeat = accumulate_trunc(sequence, length);
        sequences[length - 1].sequence = malloc(length * sizeof(char));
        memcpy(sequences[length - 1].sequence, sequence, length * sizeof(char));
        reference->length = length;
        reference->sequences = sequences;
        write_ref(*reference);
    }
}

// Compare a sequence to a reference collection and print the result.
void print_ref(char *sequence, int length, int repeat, FILE *file, reference_t reference)
{
    if (length <= reference.length)
    {
        sequence_t ref = reference.sequences[length - 1];
        if (repeat < ref.repeat)
        {
            fprintf(file, "!< %d", ref.repeat);
        }
        else if (repeat == ref.repeat)
        {
            int comp = compare(sequence, ref.sequence, length);
            if (comp < 0)
            {
                fprintf(file, "!< ");
                print(ref.sequence, length, file);
            }
            else if (comp == 0)
            {
                fprintf(file, "=");
            }
            else
            {
                fprintf(file, ">");
            }
        }
        else
        {
            fprintf(file, "> %d", ref.repeat);
        }
    }
    else
    {
        fprintf(file, "?");
        return;
    }
}

// Shared iteration state.
typedef struct
{
    int repeat;           // lowest repeat so far
    int length;           // target sequence length
    char *sequence;       // sequence buffer
    char found[MAX];      // first occurrence of lowest repeat sequence
    long long count;      // number of occurrences of lowest repeat sequence
    float last;           // last progress value
    int depth;            // current progress depth
    clock_t time;         // last progess time
    clock_t check;        // last check time
    unsigned int counter; // spinner counter
    reference_t read;     // reference sequences
    reference_t *write;   // reference sequences
} state_t;

// Initialize the iteration state.
void init(state_t *state, char *sequence, int length, reference_t read, reference_t *write)
{
    state->repeat = 0x1000000;
    state->sequence = sequence;
    state->length = length;
    state->count = 0;
    state->last = 0;
    state->depth = length - 1;
    state->time = clock();
    state->check = clock();
    state->counter = 0;
    state->read = read;
    state->write = write;
}

// Print the current progress.
void tick(state_t *state, float progress)
{
    clock_t now = clock();
    if (now > state->time + CLOCKS_PER_SEC)
    {
        int eta = (now - state->time) * (1 - progress) / (progress - state->last) / CLOCKS_PER_SEC;
        state->time = now;
        state->last = progress;
        fprintf(stderr, " %d %.2f%% %02d:%02d:%02d %dx%lld ", state->length, 100.0 * progress, eta / 3600, (eta / 60) % 60, eta % 60, state->repeat, state->count);
        print(state->found, state->length, stderr);
        fprintf(stderr, " ");
        print_ref(state->found, state->length, state->repeat, stderr, state->read);
        fprintf(stderr, "    \r");
    }
    fprintf(stderr, "%c\r", "|/-\\"[state->counter++ % 4]);
    if (now > state->check + CLOCKS_PER_SEC / 4)
    {
        state->depth += 2;
    }
    else if (now < state->check + CLOCKS_PER_SEC / 16)
    {
        state->depth--;
    }
    state->check = now;
}

// Record the current sequence if it is the lowest repeat so far.
void record(state_t *state, int length, int repeat)
{
    if (state->repeat > repeat)
    {
        state->repeat = repeat;
        state->count = 0;
        memcpy(state->found, state->sequence, sizeof(char) * length);
        merge_ref(state->sequence, length, state->write);
    }
    if (state->repeat == repeat)
    {
        state->count++;
    }
}

#define CHARSET_LENGTH (int)(sizeof(charset) / sizeof(charset[0]) - 1)

// Iterate over all sequences of the given length.
void iterate(state_t *state, int offset, int repeat, float progress, float step)
{
    for (int value = 0; value < CHARSET_LENGTH; value++)
    {
        if (offset <= state->depth)
        {
            if (value > 0)
            {
                progress += step * value;
            }
            if (offset == state->depth)
            {
                tick(state, progress);
            }
        }
        state->sequence[offset] = value;
        int next = repeat + count(state->sequence, offset);
        if (offset >= state->length - 1)
        {
            record(state, offset + 1, next);
            continue;
        }
        if (next >= state->repeat)
        {
            continue;
        }
#ifdef BOUND
        if (next > (offset + 1 - state->read.offset) * (offset - state->read.offset) / state->read.divide + state->read.add)
        {
            continue;
        }
#endif
        iterate(state, offset + 1, next, progress, step / 2);
    }
}

// Run the iteration for the given sequence length.
void run(int offset, int repeat, char *sequence, int length, reference_t read, reference_t *write)
{
    state_t state = {0};
    init(&state, sequence, length, read, write);
    iterate(&state, offset, repeat, .0, .5);
    memcpy(sequence, state.found, sizeof(char) * length);
    print(state.found, length, stdout);
    repeat = accumulate_trunc(sequence, length);
    printf(" len=%d rep=%d cnt=%lld ", length, repeat, state.count);
    print_ref(sequence, length, repeat, stdout, read);
    printf("                \n");
    fflush(stdout);
}

// Iterate over all single sequences.
void single(int length, reference_t read, reference_t *write)
{
    char sequence[MAX] = {0};
    run(1, 0, sequence, length, read, write);
}

// Base sequence used as prefix for generating sequences.
typedef struct
{
    int length;
    const char **prefixes;
} base_t;

// Iterate over all sequences with the given base sequence.
void scan(int start, base_t base, reference_t read, reference_t *write)
{
    char sequence[MAX] = {0};
    for (int length = start; length <= MAX; length++)
    {
        for (int b = 0; b < base.length; b++)
        {
            int offset = parse(base.prefixes[b], sequence, length);
            int repeat = accumulate_trunc(sequence, offset);
            run(offset, repeat, sequence, length, read, write);
        }
    }
}

// Iterate over all sequences with the given base sequence and a single character appended.
void snake(int start, int ahead, const char *base, reference_t read, reference_t *write)
{
    if (start - ahead > (int)strlen(base))
    {
        printf("!! base too short\n");
        return;
    }
    char sequence[MAX] = {0};
    parse(base, sequence, start - ahead);
    int repeat = accumulate_trunc(sequence, start - ahead);
    for (int length = start; length <= MAX; length++)
    {
        int offset = length - ahead;
        run(offset, repeat, sequence, length, read, write);
        parse(base, sequence, offset + 1);
        repeat += count(sequence, offset);
    }
}

// Check if the reference has correct repeat counts.
void check(reference_t read)
{
    for (int length = 1; length <= read.length; length++)
    {
        int repeat = accumulate_trunc(read.sequences[length - 1].sequence, length);
        if (repeat != read.sequences[length - 1].repeat)
        {
            printf("!! len=%d rep=%d ref=%d\n", length, repeat, read.sequences[length - 1].repeat);
        }
    }
}

#define LOG(length) (log(length) / log(CHARSET_LENGTH))

// Estimate the repeat count for an theoreticaly optimal sequence of the given length, based on truncated sum metric.
int optimal_sum(int length)
{
    return (double)length * (length - 1) / 2 * (length + CHARSET_LENGTH - 4) / length / (CHARSET_LENGTH - 1) - LOG(length) * (length - LOG(length) + CHARSET_LENGTH - 3) / 2;
}

// Provide a few metrics to estimate bound values.
void bound(reference_t reference)
{
    int max[MAX] = {0};
    int max_add;
    for (int length = 1; length <= reference.length; length++)
    {
        int sum = 0;
        max_add = 0;
        for (int i = 0; i < length; i++)
        {
            sum += count(reference.sequences[length - 1].sequence, i);
            if (max[i] < sum)
            {
                max[i] = sum;
            }
            int add = ceil(sum - (i + 1 - reference.offset) * (i - reference.offset) / reference.divide);
            if (max_add < add)
            {
                max_add = add;
            }
        }
        printf("%d max_add=%d\n", length, max_add);
        if (sum != reference.sequences[length - 1].repeat)
        {
            printf("!! %d: %d != %d (%d)\n", length, sum, reference.sequences[length - 1].repeat, optimal_sum(length));
        }
    }
    max_add = 0;
    for (int length = 1; length <= reference.length; length++)
    {
        int bound4 = max[length - 1] - (length - reference.offset) * (length - reference.offset - 1) / reference.divide;
        if (max_add < bound4)
        {
            max_add = bound4;
        }
        printf("%d\t%d\t%d\t%d\n", length, max[length - 1], reference.sequences[length - 1].repeat, optimal_sum(length));
    }
    printf("max_add=%d\n", max_add);
}

// Estimate the repeat count for a random sequence of the given length, using sum metric.
int estimate_sum(int length)
{
    return (length - 1) * (length + CHARSET_LENGTH - 4) / (CHARSET_LENGTH - 1) / 2;
}

// Estimate the repeat count for a random sequence of the given length, using minus one metric.
int estimate_minus(int length)
{
    return (length - 2) * (length - 3) / 4;
}

// Estimate the repeat count for a random sequence of the given length, using max metric.
int estimate_max(int length)
{
    return length * log(length * 0.3) * 1.56;
}

#define estimate(...) CONCAT(estimate_, METRIC)(__VA_ARGS__)
#define SAMPLE_COUNT 100000
// Estimate the repeat count for a random sequence of the given length.

void sample()
{
    char sequence[MAX] = {0};
    for (int length = 1; length < MAX; length++)
    {
        long long sum = 0;
        for (int s = 0; s < SAMPLE_COUNT; s++)
        {
            for (int i = 0; i < length; i++)
            {
                sequence[i] = charset[rand() / (RAND_MAX / CHARSET_LENGTH)];
                sum += count(sequence, i);
            }
        }
        int repeat = sum / SAMPLE_COUNT;
        int est = estimate(length);
        printf("len=%d rep=%d est=%d delta=%d\n", length, repeat, est, est - repeat);
    }
}

// De Bruijn recursion.
char *db(int t, int p, int n, char *a, char *s)
{
    if (t > n)
    {
        if (n % p == 0)
        {
            memcpy(s, a + 1, p * sizeof(char));
            s += p;
        }
    }
    else
    {
        a[t] = a[t - p];
        s = db(t + 1, p, n, a, s);
        for (int j = a[t - p] + 1; j < CHARSET_LENGTH; j++)
        {
            a[t] = j;
            s = db(t + 1, t, n, a, s);
        }
    }
    return s;
}

// Generate a De Bruijn sequences.
void debruijn(int ahead, reference_t read, reference_t *write)
{
    char a[MAX] = {0};
    char sequence[MAX] = {0};
    int max = LOG(MAX);
    for (int n = 1; n <= max; n++)
    {
        memset(a, 0, CHARSET_LENGTH * n * sizeof(char));
        int length = (int)pow(CHARSET_LENGTH, n);
        db(1, 1, n, a, sequence);
        print(sequence, length, stdout);
        int repeat = accumulate_trunc(sequence, length);
        printf(" len=%d rep=%d ", length, repeat);
        print_ref(sequence, length, repeat, stdout, read);
        printf("\n");
        merge_ref(sequence, length, write);
        if (ahead > 0)
        {
            int offset = length > ahead ? length - ahead : 0;
            repeat = accumulate_trunc(sequence, offset);
            run(offset, repeat, sequence, length, read, write);
        }
    }
}

// Same as the sum metric, but using wrap-around instead of truncation
long long accumulate_wrap(char *sequence, int length)
{
    long long count = 0;
    for (int k = 1; k < length; k++)
    {
        for (int offset = k - 1; offset >= 0; offset--)
        {
            for (int i = offset, j = k, l = length; l > 0 && sequence[(i + length) % length] == sequence[(j + length) % length]; i--, j--, l--)
            {
                count++;
            }
        }
    }
    return count;
}

#define EDGE trunc // trunc, wrap
#define accumulate(...) CONCAT(accumulate_, EDGE)(__VA_ARGS__)

#define MAX2 0x10000

// Provides an overview of the closed form equations for repeat counts.
void summary()
{
    char a[MAX2] = {0};
    char sequence[MAX2] = {0};
    char buffer[MAX2 * 8] = {0};
    int max = LOG(MAX2);
    for (int n = 1; n <= max; n++)
    {
        int length = (int)pow(CHARSET_LENGTH, n);
        double base = (double)length * (length - 1) / 2;
        // zeros (max)
        memset(sequence, 0, length * sizeof(char));
        long long zeros = accumulate(sequence, length);
        double zeros_closed = base * (length + 1) / 3; // trunc
        // double zeros_closed = base * length; // wrap
        // cyclic
        for (int i = 0; i < length; i++)
        {
            sequence[i] = charset[i % CHARSET_LENGTH];
        }
        long long cyclic = accumulate(sequence, length);
        double cyclic_closed = (double)length * (length - CHARSET_LENGTH) * (length * 2 - CHARSET_LENGTH + 3) / CHARSET_LENGTH / 12; // trunc
        // double cyclic_closed = zeros_closed - (double)length * length * length * (CHARSET_LENGTH - 1) / CHARSET_LENGTH / 2; // wrap
        // text
        FILE *file = fopen("file.txt", "r");
        fread(buffer, sizeof(char), (length + 7) / 8, file);
        fclose(file);
        memset(sequence, 0, length * sizeof(char));
        for (int i = 0; i < length; i++)
        {
            sequence[i] = buffer[i / 8] >> (7 - i % 8) & 1;
        }
        long long text_repeat = accumulate(sequence, length);
        // random
        long long sum = 0;
        srand(time(NULL));
        for (int s = 0; s < SAMPLE_COUNT; s++)
        {
            for (int i = 0; i < length; i++)
            {
                sequence[i] = charset[rand() / (RAND_MAX / CHARSET_LENGTH)];
            }
            sum += accumulate(sequence, length);
        }
        long long rnd = sum / SAMPLE_COUNT;
        double rnd_closed = base * (length + CHARSET_LENGTH - 4) / length / (CHARSET_LENGTH - 1); // trunc
        // double rnd_closed = base / (CHARSET_LENGTH - 1); // wrap
        // de Bruijn (min)
        memset(a, 0, CHARSET_LENGTH * n * sizeof(char));
        db(1, 1, n, a, sequence);
        long long repeat = accumulate(sequence, length);
        double repeat_closed = rnd_closed - LOG(length) * (length - LOG(length) + CHARSET_LENGTH - 3) / 2; // trunc
        // double repeat_closed = rnd_closed - (LOG(length) * length / 2); // wrap
        // summary
        printf("len=%d min=%lld,%.2f rnd=%lld,%.2f text=%lld cyc=%lld,%.2f max=%lld,%.2f\n", length, repeat, repeat_closed, rnd, rnd_closed, text_repeat, cyclic, cyclic_closed, zeros, zeros_closed);
        // printf("len=%d min=%.4f->%.2f rnd=%.4f~=%.2f text=%.4f cyc=%.4f->%.2f max=%.4f==1\n", length, repeat/base, 1./(CHARSET_LENGTH - 1), rnd/base, 1./(CHARSET_LENGTH - 1), text_repeat/base, cyclic/zeros_closed, 1 - (CHARSET_LENGTH - 1.0) / CHARSET_LENGTH, zeros/zeros_closed);
    }
}

// Computes the repeat count for a binary (or text) sequence read from a file.
void file(const char *filename, int max)
{
    FILE *file = fopen(filename, "r");
    fseek(file, 0L, SEEK_END);
    int size = ftell(file) + 1;
    size = size < (max + 7) / 8 ? size : (max + 7) / 8;
    char *buffer = malloc(size * sizeof(char));
    fseek(file, 0L, SEEK_SET);
    size = fread(buffer, sizeof(char), size, file);
    fclose(file);
    int length = max < size * 8 ? max : size * 8;
    char *sequence = malloc(length * sizeof(char));
    for (int i = 0; i < length; i++)
    {
        sequence[i] = buffer[i / 8] >> (7 - i % 8) & 1;
    }
    free(buffer);
    long long repeat = accumulate(sequence, length);
    free(sequence);
    long long random = (long long)length * (length - 1) / 2;
    double relative = (double)repeat / random;
    printf("len=%d rep=%lld rel=%.3f\n", length, repeat, relative);
}

int main()
{
    reference_t read = {0, 1., 1., 1., NULL};
    reference_t write = {0, 1., 1., 1., NULL};
    read_ref(&read);
    read_ref(&write);

    single(30, read, &write);
    // scan(1, (base_t){ 1, (const char *[1]){ "0" } }, read, &write);
    // scan(42, (base_t){ 6, (const char *[6]){ "00000010", "00000100", "00001000", "00010000", "00100000", "01000000" } }, read, &write);
    // snake(42, 30, "000000010000", read, &write);
    // check(read);
    // bound(read);
    // sample();
    // debruijn(30, read, &write);
    // summary();
    // file("file.txt", 0x800);

    return 0;
}
