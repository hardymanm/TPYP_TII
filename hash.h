#define HASH_TABLE_MAX_SIZE 100000000
typedef struct HashNode_Struct HashNode;

struct HashNode_Struct
{
    char* sKey;
    double nValue;
    HashNode* pNext;
};

void hash_table_init();
unsigned int hash_table_hash_str(const char* skey);
void hash_table_insert(const char* skey, double nvalue);
void hash_table_remove(const char* skey);
HashNode* hash_table_lookup(const char* skey);
void hash_table_print();
void hash_table_release();
